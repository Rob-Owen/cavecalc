"""Performs geochemical calculations and handles IO with IPhreeqc.

This module contains the core functionality of the Cavecalc model. All
geochemical calculations are defined here. Some are also calculated here - 
others are delegated to IPhreeqc.

Contains classes:
    Solution    - solution chemistry declaration and reactions
    Gas         - gas phase declaration and calculations
    Carbonate   - carbonate phase definition and calculations
    Simulator   - handles IO with IPhreeqc
"""

import sys
import os
import math
import re
import numpy as np
from decimal import *
import cavecalc.util as ccu
from cavecalc.data.phreeqc_templates import *
from cavecalc.setter import SettingsObject
from copy import copy, deepcopy
import inflection

## phreeqpy mode (phreeqpy_dll, API_com or phreeqpy_com):
# mode may affect speed - author found phreeqpy_dll fastest but mileage may
# vary depending on system. 
#   - API_com requires Windows and a full IPhreeqcCOM install
#   - phreeqpy_dll requires phreeqpy with a suitable IPhreeqc .dll/.so file
#   - phreeqpy_com requires Windows, a full IPhreeqcCOM install and phreeqpy
if sys.platform == 'win32':
    MODE = 'API_com'
else:
    MODE = 'phreeqpy_dll'

# edit MODE here to force behaviour

if MODE == 'phreeqpy_com':
    import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
elif MODE == 'phreeqpy_dll':
    import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
elif MODE == 'API_com':
    from win32com.client import Dispatch
else:
    raise ValueError('Mode "%s" is not defined.' % MODE)
    
# Global parameters
EQ_STEP = True                      #   run initial isotope equilibration step
PHREEQC_TOLERANCE = 0.1             #   tolerance on phreeqc percent error
PRESSURE = 1.001                    #   atmospheric pressure (bar)
BEDROCK_PHASE_QZ_LEVEL = 1e-10      #   precision of bedrock phase declaration

# Class definitions
class Solution(object):
    """Performs calculations for aqueous reactions.
    
    Solution objects handle aqeuous chemistry reactions. Most methods return
    PHREEQC input strings that may be passed to IPhreeqc via a Simulator
    object. The Solution object does not store any data concerning the present
    solution state - that is handled by a Simulator object.
    
    This is usually instantiated by a Simulator object.
    
    Selected Methods (return reaction strings):
        Solution.set_soil_chem
        Solution.open_diss_rxn
        Solution.open_diss_rxn_with_reprecip
        Solution.wri_rxn
        Solution.full_degas_rxn
        Solution.degas_rxn
    """
    
    def __init__(self, Sim):
        """Initialise a solution object.
        
        Args:
            Sim: The Simulator object to provide input settings and handle 
            IO with IPhreeqc.
        """

        self.s = Sim
        self.G = Gas(Sim=self.s)
        self.init_O = None # initial solution O concentration
        self.init_C = None # initial solution C concentration
        self.init_pH = None # initial solution pH (improved estimate)
        self.init_13c_aq = None # DIC d13c_aq
        
        self._get_init_co()
        self._get_init_d13c()
     
    @property
    def calcite_string(self):
        """Returns an empty calcite solid solution string for Iphreeqc."""
        
        if hasattr(self, '_calcite_string'):
            return "\nUSE SOLID_SOLUTION 1\n"
        else:
            # create solid solution template
            self._calcite_string = CALCITE_SS
            if self.s.settings['soil_Mg'] > 0 or \
               self.s.settings['bedrock_MgCa'] > 0 or \
               self.s.settings['bedrock_mineral'].capitalize() == 'Dolomite':
                self._calcite_string += SS_MG   # add Mg-bearing members
            if self.s.settings['soil_Sr'] > 0 or \
               self.s.settings['bedrock_SrCa'] > 0:
                self._calcite_string += SS_SR   # add Sr-bearing members                
            if self.s.settings['soil_Ba'] > 0 or \
               self.s.settings['bedrock_BaCa'] > 0:
                self._calcite_string += SS_BA   # add Ba-bearing members
            return self._calcite_string

    @property
    def diss_gas_string(self):
        """Returns a IPhreeqc GAS_PHASE definition for the soil gas specified
        by user inputs. This is used for bedrock dissolution reactions.
        """
    
        if hasattr(self, '_diss_gas'):
            return self._diss_gas.get_gasphase_string()
            
        DissGas = Gas(Sim=self.s)
        DissGas.set_initial_gas()
        self._diss_gas = DissGas
        return self._diss_gas.get_gasphase_string()
            
    @property
    def degas_string(self):
        """Prepare PHREEQC input for degassing a small amount of CO2.
        
        Degassing is unidirectional and uses the REACTION keyword in its
        output. Isotope fractionation factors are taken from the database file.
        The fraction of aqueous CO2 removed is defined in the model settings.
        
        Returns:
            (pq_str, desc_str) for adding to the iphreeqc buffer.
        """
        
        
        co2_change = ( self.s.settings['co2_decrement'] *
                        self.s.get('m_CO2(mol/kgw)') *
                        self.s.get('mass_H2O') )
        
        a = self.G.calculate_eq_co2(co2_change)
        
        fmt = { 'molefrac_co2'  :   a[0],
                'molefrac_13C'  :   a[1],
                'molefrac_18O'  :   a[3],
                'molefrac_14C'  :   a[2],
                'gas_moles'     :   -1*co2_change}
                
        #degas_input = CO2_OUT.format(**fmt) # include radiocarbon (unnecessary)
        degas_input = CO2_OUT2.format(**fmt) # no radiocarbon (greater stability)
        return degas_input
        
    def _get_init_co(self):
        """Calculate initial solution pH, C and O concentrations.

        This sets up a separate Simulator object to to determine dissolved
        C and O concentrations in the initial solution. This step is
        performed before self._set_init_chem() as PHREEQC cannot determine
        gas-solution equilibrium at the same time as declaring solution isotope
        compositions.
        """
    
        # get settings
        fmt = { 'temp'  :   self.s.settings['temperature'],
                'ph'    :   self.s.settings['soil_pH'],
                'ca'    :   self.s.settings['soil_Ca'],
                'mg'    :   self.s.settings['soil_Mg'],
                'sr'    :   self.s.settings['soil_Sr'],
                'ba'    :   self.s.settings['soil_Ba'],
                'cl'    :   2*( self.s.settings['soil_Ca'] + 
                                self.s.settings['soil_Mg'] +
                                self.s.settings['soil_Sr'] +
                                self.s.settings['soil_Ba'] ),
                'co2_si':   math.log10(PRESSURE * 1e-6 * self.s.settings['init_pCO2']) }
        
        if self.s.settings['init_O2'] == 0:
            fmt['o2_si'] = -10
        else:
            fmt['o2_si'] = math.log10(PRESSURE * self.s.settings['init_O2'])
        
        # create a dummy PHREEQC solution in a new Simulator to calculate C and
        # O concentrations in initial solution
        dummy_chemistry = DUMMY_CHEM_CO.format(**fmt)
        dummy_settings = self.s.settings.copy()
        dummy_settings['phreeqc_log_file_name'] = 'log_{}_init.phr'
        dummy_simulator = Simulator(dummy_settings, self.s.id)
        dummy_simulator.ipq_buffer(dummy_chemistry, 'Calculate init C/O concs')
        dummy_simulator.ipq_exec()
        
        if self.s.settings['init_O2'] == 0:
            self.init_O = 0
        else:
            self.init_O = dummy_simulator.get('O(mol/kgw)')*1000 # initial [O]
            
        self.init_C = dummy_simulator.get('C(mol/kgw)')*1000 # initial [C]
        self.init_pH = dummy_simulator.get('pH')    # improved estimate of pH
        
    def _get_init_d13c(self):
        """Calculate the initial solution d13C.
        
        This requires it's own calculation step because init_d13C is given as
        a gas phase composition - we must calculate the initial solution ratio.
        
        To calculate the d13C of the solution based on init_d13C (a gas phase 
        composition), we do the equilibration between a large gas volume with
        init_d13C and a small amount of solution, initially also at init_d13C.
        From this the 'correct' initial d13C may be calculated.
        """
        
        t = self.s.settings['temperature']
        fmt = { 'temp'  :   t,
                'ph'    :   self.init_pH,
                'ca'    :   self.s.settings['soil_Ca'],
                'mg'    :   self.s.settings['soil_Mg'],
                'sr'    :   self.s.settings['soil_Sr'],
                'ba'    :   self.s.settings['soil_Ba'],
                'cl'    :   2*( self.s.settings['soil_Ca'] + 
                                self.s.settings['soil_Mg'] +
                                self.s.settings['soil_Sr'] +
                                self.s.settings['soil_Ba'] ),
                'c'     :   self.init_C,
                'o'     :   self.init_O,
                'd18O'  :   self.s.settings['atm_d18O'],
                'd13C'  :   self.s.settings['init_d13C'], # guesstimate
                'R14C'  :   self.s.settings['init_R14C'],
                'd44Ca' :   self.s.settings['soil_d44Ca']}
        
        dummy_chemistry = (DUMMY_CHEM_D13C.format(**fmt), 
                            'soln with guesstimate d13C')

        self.G.set_initial_gas()
        self.G.volume = 1000 # arbitrary large volume
        self.G.name = 'Equilibration Gas'
        gas = self.G.get_gasphase_string()

        dummy_gas = (gas, 'equilibrate with high vol gas /w correct d13C')
        
        dummy_settings = self.s.settings.copy()
        dummy_settings['phreeqc_log_file_name'] = 'log_{}_d13c.phr'
        dummy_simulator = Simulator(dummy_settings, self.s.id)
        
        d_conv = [self.s.settings['init_d13C']]
        ph_conv = [self.init_pH]
        c_conv = [self.init_C]
        
        dummy_simulator.ipq_buffer([dummy_chemistry, dummy_gas])
        dummy_simulator.ipq_exec()
        
        # note output d18O will be slightly incorrect. Limitation of PHREEQC.
        self.init_13c_aq = dummy_simulator.get('I_R(13C)')
        
        # if a solution d13C has been explicitly stated, use that value
        if self.s.settings['init_solution_d13c'] != 'n/a':
            self.init_13c_aq = self.s.settings['init_solution_d13c']
    
    def set_soil_chem(self):
        """Add Initial Solution Chemistry to the phreeqc buffer. Initial 
        Solution Chemistry is taken from model input parameters.
        
        This method normally requires _get_init_co() and _get_init_d13C to be
        run first to set the initial solution C, O and d13C values.
        
        Returns:
            (pq_str, desc_str) for adding to the iphreeqc buffer.
        """
        
        # get initial chemistry parameters (assume Cl- for charge balance)
        init_r14c = ccu.pmc_denormalise(   self.s.settings['init_R14C'], 
                                           self.init_13c_aq    )
        
        fmt = { 'temp'  :   self.s.settings['temperature'],
                'ph'    :   self.init_pH,
                'ca'    :   self.s.settings['soil_Ca'],
                'mg'    :   self.s.settings['soil_Mg'],
                'sr'    :   self.s.settings['soil_Sr'],
                'ba'    :   self.s.settings['soil_Ba'],
                'cl'    :   2*( self.s.settings['soil_Ca'] + 
                                self.s.settings['soil_Mg'] +
                                self.s.settings['soil_Sr'] +
                                self.s.settings['soil_Ba'] ),
                'c'     :   self.init_C,
                'o'     :   self.init_O,
                'd18O'  :   self.s.settings['atm_d18O'],
                'd13C'  :   self.init_13c_aq,
                'R14C'  :   init_r14c,
                'd44Ca' :   self.s.settings['soil_d44Ca']}
        
        self.initial_chemistry = INITIAL_CHEM.format(**fmt)
                                     
        if EQ_STEP:
            return self.initial_chemistry + WATER_EQUILIBRATE, 'initial_water'
        else:
            return self.initial_chemistry, 'initial_water'   
            
    def open_diss_rxn(self, Carbonate):
        """Prepare PHREEQC input for open-system dissolution.
        
        Details of the bedrock and gas phase present are extracted from 
        self.sim.settings. By varying gas_volume from 0 (closed system) up to 
        <very big numbers>, open vs closed system dissolution may be 
        quantified.
        
        open_diss_rxn with 0 gas volume is equivalent to closed system
        dissolution.
        
        Args:
            Carbonate: A Carbonate object defining the solid phase.
        Returns:
            (pq_str, desc_str) for adding to the iphreeqc buffer.
        """        
        
        rxn_string = ""
        
        if not self.s._bedrock_set:
            rxn_string += Carbonate.phase_definition + '\n'
            self.s._bedrock_set = True
        
        rxn_string += Carbonate.bedrock_equilibrium1()
        rxn_string += '\n' + self.diss_gas_string
        return rxn_string, 'dissolve bedrock'
                    
    def open_diss_rxn_with_reprecip(self, Carbonate):
        """Prepare PHREEQC input for dissolution with reprecipitation.
        
        Identical to open_diss_rxn but allows coeval precipitation of calcite
        while bedrock is dissolving. This simulates incongruent bedrock 
        dissolution. 
        
        Note that this allows extreme amounts of reprecipitation. For more
        controlled re-precipitation, use open_diss_rxn followed by wri (water-
        rock interactions).

        Args:
            Carbonate: A Carbonate object defining the solid phase.
        Returns:
            (pq_str, desc_str) for adding to the iphreeqc buffer.
        """

        rxn_string = ""
        
        if not self.s._bedrock_set:
            rxn_string += Carbonate.phase_definition + '\n'
            self.s._bedrock_set = True
        
        rxn_string += Carbonate.bedrock_equilibrium2()
        rxn_string += '\n' + self.diss_gas_string
        rxn_string += '\n' + self.calcite_string
        return rxn_string, 'dissolve bedrock'
    
    def wri_rxn(self, Carbonate, m_bedrock=0, m_pyrite=0):
        """Prepare phreeqc input for water-rock interactions.
        
        The resulting Iphreeqc command allows the current solution to dissolve
        a set amount of bedrock. No capacity is included for gas exchange or
        calcite precipitation - these reactions should be added separately.
        
        This differs from open_diss_rxn_with_reprecip, which allows
        re-precipitating in the presence of gas.  

        Args:
            Carbonate: The Carbonate object to derive bedrock definitions from.
            m_bedrock: Moles of bedrock present. Default is 0.
            m_pyrite: Moles of pyrite present. Default is 0.
        Returns:
            (pq_str, desc_str) for adding to the iphreeqc buffer.
        """
        
        rxn_string = ""
        
        if not self.s._bedrock_set:
            rxn_string += Carbonate.phase_definition + '\n'
            self.s._bedrock_set = True
        
        rxn_string += Carbonate.bedrock_equilibrium2( m_bedrock=m_bedrock, 
                                                      m_pyrite=m_pyrite)
        # rxn_string += '\n' + self.calcite_string
        # return rxn_string, 'water-rock equilibration & precipitation'
        return rxn_string, 'water-rock equilibration'
        
    def eq_degas_rxn(self, EqGas=None):
        """Prepare PHREEQC input for solution-gas CO2 equilibration.

        Args:
            EqGas: A Gas object to equilibrate the solution with. If none is
                   provided, cave air gas is used.
        Returns:
            (pq_str, desc_str) for adding to the iphreeqc buffer.
        """
        
        if Gas is None:
            EqGas = Gas(Sim=self.s)
            EqGas.set_cave_air()
            
        rxn_string = EqGas.get_gasphase_string()
        return rxn_string, 'Equilibrium Degassing'
    
    def kinetic_degas_rxn(self):
        """Prepare PHREEQC input for degassing a small amount of CO2.
        
        Degassing is unidirectional and uses the REACTION keyword in its
        output. Isotope fractionation factors are taken from the database file.
        The fraction of aqueous CO2 removed is defined in the model settings.
        
        Returns:
            (pq_str, desc_str) for adding to the iphreeqc buffer.
        """
        
        return self.degas_string, 'Kinetic Degassing'
        
    def will_it_degas(self, pCO2, printout=False):
        """Test whether the solution will degas at a given pCO2.
        
        Args:
            pCO2: gas CO2 concentration (ppmv)
            printout (logical): If true, prints the value of Solution pCO2
        Returns:
            True/False
        """
        
        temp         = self.s.get('temp(C)') + 273.15
        water        = self.s.get('mass_H2O')
        soln_CO2     = self.s.get('m_CO2(mol/kgw)')
        log_k        = self.s.reader.get_k_values([r'CO2\(g\)'], temp)['analytic_value']
        eq_CO2       = soln_CO2 / 10**log_k
        gas_CO2      = pCO2 * 1e-6   
        
        self.co2_change = self.s.settings['co2_decrement'] * soln_CO2 * water
        
        if printout:
            print("Water pCO2 (ppmv):\t%f"    %(eq_CO2 * 1e6))
            print("Gas phase pCO2 (ppmv):\t%f" % pCO2)
        
        if eq_CO2 < gas_CO2:
            if __name__ == '__main__':
                print("Water pCO2 (mol/kgw):\t%f" % soln_CO2)
                print("Water pCO2 (ppmv):\t%f"    %(eq_CO2 * 1e6))
                print("Gas phase pCO2 (ppmv):\t%f" % pCO2)
                print("Degassing complete.")
            return False
        else:
            return True
    
    def precipitate_rxn(self):
        """Prepare PHREEQC input for calcite precipitation.
        
        Precipitation will occur until the solution reaches saturation with
        respect to calcite. Radiocarbon fractionation is not included as the
        very small amounts present may cause numerical instability.
        
  
        Returns:
            (pq_str, desc_str) for adding to the iphreeqc buffer.
        """

        return self.calcite_string, "CaCO3_precipitation"
    
class Carbonate(object):
    """Defines bedrock chemistry.
    
    Solid phase isotope assemblages are not easily expressed natively in
    PHREEQC. This class allows isotopically-distinct calcite and dolomite
    mineralogies to be defined as PHREEQC phases for the purposes of
    dissolution reactions.
    
    Carbonate objects take their geochemical definitions from the settings
    of a Simulator object.
    
    Selected Attributes:
        phase_definition - bedrock phase definition string

    Selected Methods:
        bedrock_equilibrium1 - bedrock dissolution string
        bedrock_equilibrium2 - bedrock dissolution string (alternative)
    """
    
    def __init__(self, Sim):
        """Constructs PHREEQC input for isotope-defined bedrock.
        
        Args:       
            Sim: The Simulator object to draw input from.
        
        Returns:
            Nothing
        """
        
        self.s = Sim
        
        if self.s.settings['bedrock_mineral'].capitalize() == 'Calcite':
            # self._define_bedrock('Calcite')
            self._get_calcite_phase()
        elif self.s.settings['bedrock_mineral'].capitalize() == 'Dolomite':
            self._get_dolomite_phase()
        else:
            raise TypeError("Bedrock mineralogy not recognised: %s" %
                                self.s.settings['bedrock_mineral'])
    
        # self.check_charge_balance(self.diss_reaction)
       
    @property
    def phase_definition(self):
        if hasattr(self, '_pd'):
            return self._pd
            
        self._pd = BEDROCK_PHASES.format( self.pq_phase_entry )
        return self._pd
            
    def bedrock_equilibrium1(self, m_bedrock=None, m_pyrite=None):
        """Returns a PHREEQC input string for bedrock dissolution.

        This method should not be used for incongruent dissolution reactions as
        the name 'Calcite' may not be used twice in the same PHREEQC step. In
        such cases use bedrock_equilibrium2.
        
        Args:
            m_bedrock: Moles of bedrock present during dissolution. Default
                       uses model input parameters.
            m_pyrite: Moles of pyrite present during dissolution. Default
                      uses model input parameters.
        Returns:
            A PHREEQC input string specifying a bedrock dissolution reaction.
        """
        
        if m_bedrock is None:
            m_bedrock = self.s.settings['bedrock']
        if m_pyrite is None:
            m_pyrite = self.s.settings['bedrock_pyrite']
              
        eq_string = BEDROCK_DISS1.format( moles_bedrock = m_bedrock,
                                          moles_pyrite = m_pyrite ) 
        return eq_string

    def bedrock_equilibrium2(self, m_bedrock=None, m_pyrite=None):
        """Returns a PHREEQC input string for bedrock dissolution.
        
        Alternative to bedrock_equilibrium1 when Calcite is present during the
        same reaction step (most commonly during incongruent dissolution
        calculations).
        
        Args:
            m_bedrock: Moles of bedrock present during dissolution. Default
                       uses model input parameters.
            m_pyrite: Moles of pyrite present during dissolution. Default
                      uses model input parameters.
        Returns:
            A PHREEQC input string specifying a bedrock dissolution reaction.
        """
            
        if m_bedrock is None:
            m_bedrock = self.s.settings['bedrock']
        if m_pyrite is None:
            m_pyrite = self.s.settings['bedrock_pyrite']
              
        eq_string = BEDROCK_DISS2.format( moles_bedrock = m_bedrock,
                                          moles_pyrite = m_pyrite ) 
        return eq_string

    def _get_calcite_phase(self):
        
        def dcqz(decimal):
            d2 = decimal.quantize( Decimal(str(BEDROCK_PHASE_QZ_LEVEL)),
                                   rounding=ROUND_HALF_UP)
            return d2
            
        bed_44Ca = self.s.settings['bedrock_d44Ca']
        bed_13C = self.s.settings['bedrock_d13C']
        bed_18O = self.s.settings['bedrock_d18O']
        MgCa = Decimal(self.s.settings['bedrock_MgCa'] * 0.001)
        SrCa = Decimal(self.s.settings['bedrock_SrCa'] * 0.001)
        BaCa = Decimal(self.s.settings['bedrock_BaCa'] * 0.001)
        a = Decimal(str(self.s.stnd44Ca * (1 + 0.001*bed_44Ca))) # 44/40 ratio
        b = Decimal(str(self.s.stnd13C * (1 + 0.001*bed_13C))) # 13/12 ratio
        c = Decimal(str(self.s.stnd18O * (1 + 0.001*bed_18O))) # 18/16 ratio
    
        # metal totals
        tCa = 1 / (1+MgCa+SrCa+BaCa)
        tMg = dcqz(MgCa * tCa)
        tSr = dcqz(SrCa * tCa)
        tBa = dcqz(BaCa * tCa)
        
        # isotope totals
        Ca44 = dcqz(a/(1+a) * tCa) 
        Ca40 = dcqz(tCa - Ca44)
        C13 = dcqz(b/(1+b))
        C12 = dcqz(1 - C13)
        O18 = dcqz(c/(1+c) * 3)
        O16 = dcqz(3 - O18)
        
        # construct phase definition
        metals = "Ca{Ca40}[44Ca]{Ca44}".format(Ca40=Ca40, Ca44=Ca44)
        d = "{}Ca+2 + {}[44Ca]+2 + {}[13C]O3-2 + {}CO2[18O]-2 + {}CO3-2"
        d = d.format(Ca40, Ca44, C13, O18, C12-O18)
        
        if MgCa > 0:
            metals += "Mg{Mg}".format(Mg=tMg)
            d += " + {}Mg+2".format(tMg)
        if SrCa > 0:
            metals += "Sr{Sr}".format(Sr=tSr)
            d += " + {}Sr+2".format(tSr)
        if BaCa > 0:
            metals += "Ba{Ba}".format(Ba=tBa)
            d += " + {}Ba+2".format(tBa)
            
        carbonate = "C{}[13C]{}O{}[18O]{}".format( C12, C13, O16, O18 )
        s = metals + carbonate
         
        v = self.s.reader.get_k_values([r'CaCO3'])
        
        "{}Mg+2 + {}Sr+2 + {}Ba+2"
        self.diss_reaction = s + ' = ' + d
        self.pq_phase_entry = self.diss_reaction + \
                              '\n\tlog_k\t' + str(v['log_k']) +  \
                              '\n\tdelta_h\t' + str(v['delta_h']) + \
                              '\n\t' + str(v['analytic_line'])
    
    def _get_dolomite_phase(self):
        
        def dcqz(decimal):
            d2 = decimal.quantize( Decimal(str(BEDROCK_PHASE_QZ_LEVEL)),
                                   rounding=ROUND_HALF_UP)
            return d2

        bed_44Ca = self.s.settings['bedrock_d44Ca']
        bed_13C = self.s.settings['bedrock_d13C']
        bed_18O = self.s.settings['bedrock_d18O']
        SrCa = Decimal(self.s.settings['bedrock_SrCa'] * 0.001)
        BaCa = Decimal(self.s.settings['bedrock_BaCa'] * 0.001)
        a = Decimal(str(self.s.stnd44Ca * (1 + 0.001*bed_44Ca))) # 44/40 ratio
        b = Decimal(str(self.s.stnd13C * (1 + 0.001*bed_13C))) # 13/12 ratio
        c = Decimal(str(self.s.stnd18O * (1 + 0.001*bed_18O))) # 18/16 ratio
    
        # metal totals
        tCa = 2 / (2+SrCa+BaCa)
        tMg = dcqz(2 / (2+SrCa+BaCa))
        tSr = dcqz(SrCa * tCa)
        tBa = dcqz(BaCa * tCa)
        
        # isotope totals
        Ca44 = dcqz(a/(1+a) * tCa)
        Ca40 = dcqz(tCa - Ca44)
        C13 = dcqz(2*b/(1+b))
        C12 = dcqz(2 - 2*b/(1+b))
        O18 = dcqz(c/(1+c) * 6)
        O16 = dcqz(6 - c/(1+c) * 6)
        
        # construct phase definition
        metals = "Ca{Ca40}[44Ca]{Ca44}".format(Ca40=Ca40, Ca44=Ca44)
        metals += "Mg{}".format(tMg)
         
        d = "{}Ca+2 + {}[44Ca]+2 + {}Mg+2 + {}[13C]O3-2 + {}CO2[18O]-2 + {}CO3-2"
        d = d.format(Ca40, Ca44, tMg, C13, O18, C12-O18)
        
        if SrCa > 0:
            metals += "Sr{}".format(tSr)
            d += " + {}Sr+2".format(tSr)
        if BaCa > 0:
            metals += "Ba{}".format(tBa)
            d += " + {}Ba+2".format(tBa)
            
        carbonate = "C{C12}[13C]{C13}O{O16}[18O]{O18}".format(
                       C12=C12, C13=C13, O16=O16, O18=O18
                       )
             
        s = metals + carbonate

        
        v = self.s.reader.get_k_values([r"CaMg\(CO3\)2"])
        
        self.diss_reaction = s + ' = ' + d
        self.pq_phase_entry = self.diss_reaction + \
                              '\n\tlog_k\t' + str(v['log_k']) +  \
                              '\n\tdelta_h\t' + str(v['delta_h'])
        
    @classmethod
    def check_charge_balance(reaction_string):
        """Checks a PHREEQC dissolution string for charge balance.
        
        The string checked must be of the form SOLID = A + B + C where A, B, C
        are aqueous ions and SOLID is a single solid phase. This form is common
        in PHASES definitions where a mineral dissolution reaction is
        specified.
        
        This method is primarly used as a debugging tool to check the strings
        created by _get_calcite_phase() and _get_dolomite_phase() for
        isotopically and trace-element distinct bedrock phases.
        
        Args:
            reaction_string (str) A phreeqc-valid mineral dissolution reaction
                string.
        """

        # Carbonate oxidation states
        E = {   'Ca'    :   2,
                '[44Ca]':   2,
                'C'     :   4,
                '[13C]' :   4,
                'O'     :   -2,
                '[18O]' :   -2,
                'Mg'    :   2,
                'Sr'    :   2,
                'Ba'    :   2
            }
        
        print("\n" + reaction_string)
        lhs, rhs = [s.strip() for s in reaction_string.split('=')]
        
        a = []
        b = []
        
        r = re.match(r"(\[?[0-9]*[A-Z][a-z]*?\]?)([0-9]\.[0-9]*)", lhs)
        while r:
            a.append(r.groups())
            lhs = lhs.replace(r.group(0),'')
            r = re.match(r"(\[?[0-9]*[A-Z][a-z]*?\]?)([0-9]\.[0-9]*)", lhs)
        
        for species in [s.strip() for s in rhs.split(' + ')]:
            r = re.compile(r"([0-9]\.[0-9]*)(.*?)([\+\-][0-9])")
            m = re.match(r, species)
            b.append(m.groups())
        
        lhs_charge = 0
        rhs_charge = 0
        
        for m, n in a:
            lhs_charge += E[m] * Decimal(n)
        
        for m, _, n in b:
            rhs_charge += Decimal(m) * Decimal(n)
        
        print(lhs_charge)
        print(rhs_charge)
        if lhs_charge == rhs_charge:
            print("Equation is charge balanced!")
        else:
            print("Equation is NOT charge balanced :(")

class Gas(object):
    def __init__(self, pCO2=0, O2=0, d13C=0, d18O=0, R14C=100, 
                 temperature=None, volume=0, name='Gas', Sim=None):
                 
            self.pCO2 = pCO2
            self.O2 = O2
            self.d13C = d13C
            self.d18O = d18O
            self.R14C = R14C
            self.volume = volume
            self.name = name
            self.s = Sim
             
            if temperature is None:
                self.temperature = self.s.settings['temperature']
            else:
                self.temperature = temperature

    def set_soil_gas(self):
        """Set self to hold soil gas end-member chemistry."""
        
        self.pCO2 = self.s.settings['soil_pCO2']
        self.d13C = self.s.settings['soil_d13C']
        self.O2 = self.s.settings['soil_O2']
        self.d18O = self.s.settings['atm_d18O']
        self.R14C = self.s.settings['soil_R14C']
        self.volume = self.s.settings['gas_volume']
        self.name = "Soil Gas End-member"
    
    def set_atmosphere(self):
        """Set self to hold atmospheric chemistry."""
    
        self.pCO2 = self.s.settings['atm_pCO2']
        self.d13C = self.s.settings['atm_d13C']
        self.O2 = self.s.settings['atm_O2']
        self.d18O = self.s.settings['atm_d18O']
        self.R14C = self.s.settings['atm_R14C']
        self.volume = 1000 # arbitrary large volume by default
        self.name = "Atmosphere End-member"
        
    def set_initial_gas(self):
        """Set self to hold initial (mixed) soil gas chemistry."""
    
        self.pCO2 = self.s.settings['init_pCO2']
        self.d13C = self.s.settings['init_d13C']
        self.O2 = self.s.settings['init_O2']
        self.d18O = self.s.settings['atm_d18O']
        self.R14C = self.s.settings['init_R14C']
        self.volume = self.s.settings['gas_volume']
        self.name = "Mixed Initial Gas"
    
    def set_cave_air(self):
        """Set self to hold cave air chemistry."""
        
        self.pCO2 = self.s.settings['cave_pCO2']
        self.d13C = self.s.settings['cave_d13C']
        self.O2 = self.s.settings['cave_O2']
        self.d18O = self.s.settings['cave_d18O']
        self.R14C = self.s.settings['cave_R14C']
        self.volume = self.s.settings['cave_air_volume']
        self.name = "Cave Air"
             
    def calculate_eq_co2(self, co2_out):
        """Calculate CO2(g) in equilibrium with the current CO2(aq).
        
        Calculate co2 isotopologue abundances in a gas exolved from a
        solution. Equilibrium is with the CO2(aq) pool only - DIC equilibration
        is not considered. Isotope ratios are calculated via mass balance using
        fractionation factors from the database file.
        
        Args:
            co2_out: moles of CO2 to be removed from the solution.
        Returns:
            (c12o_2, c13o_2, c14o_2, c12o18o) - tuple of relative abundances
        """

        CO2_aq  = self.s.get('m_CO2(mol/kgw)')
        water   = self.s.get('mass_H2O')
        CO2_13C = self.s.get('I_R(13C)_CO2(aq)')    
        CO2_14C = self.s.get('I_R(14C)_CO2(aq)')
        CO2_18O = self.s.get('I_R(18O)_CO2(aq)')
        temp    = self.s.get('temp(C)') + 273.15
        n_co2   = CO2_aq * water
      
        assert n_co2 > co2_out
        
        # get CO2(g)-CO2(aq) fractionation factors
        frac13C = self.s.reader.get_alpha( '13C','CO2(g)/CO2(aq)', temp ) 
        frac18O = self.s.reader.get_alpha( '18O','CO2(g)/CO2(aq)', temp )
        
        # calculate normalised isotope ratios
        c13 = 1 + 0.001*CO2_13C
        o18 = 1 + 0.001*CO2_18O
        
        # calculate gas isotopes via mass balance
        c13_g = c13 * n_co2 / (co2_out + frac13C * (n_co2-co2_out))
        o18_g = o18 * n_co2 / (co2_out + frac18O * (n_co2-co2_out))
        d13C_gas = 1000*(c13_g - 1)
        d18O_gas = 1000*(o18_g - 1)
        
        # get absolute gas isotopologue abundances
        return self.co2_isotopologues(d13C=d13C_gas, d18O=d18O_gas, 
                                      R14C=CO2_14C)
                        
    def co2_isotopologues(self, d13C=None, d18O=None, R14C=None):
        """Calculate CO2 isotopologue relative abundances.
        
        Calculate relative abundances of CO2, [13C]O2, CO[18O] and [14C]O2
        based on delta/ratio values. The values sum to 1. This method uses
        the absolute isotopic composition of the standards from the phreeqc
        database and then solves the mass balance equations for each isotope.
        
        Args:
            d13C: CO2(g) d13C
            d18O: CO2(g) d18O
            R14C: CO2(g) R14C as a pmc value
        Returns:
            (c12o_2, c13o_2, c14o_2, c12o18o) - tuple of molar abundances
        """
        
        if d13C is None:
           d13C = self.d13C
        if d18O is None:
           d18O = self.d18O
        if R14C is None:
           R14C = self.R14C
        
        # calculate absolute isotope ratios (14C from Stuvier & Polach, 1977)
        s13C = self.s.stnd13C
        s18O = self.s.stnd18O
        # s14C = self.s.stnd14C
        
        r13 = s13C * (1 + 0.001 * d13C) # 13C/12C ratio
        r18 = s18O * (1 + 0.001 * d18O) # 18O/16O ratio
        r14 = ccu.pmc_2_c14(R14C, d13C, self.s.stnd14C)
        
        # r14 = s14C * 0.01 * R14C * math.pow((1+0.001*float(d13C))/0.975,2)
        # mole fraction 14C
        
        ## system of equations to solve:
        # r13C*CO2 - 13CO2 + 0*14CO2 + r13C*CO18O = 0 (13C mass balance)
        # r14C*CO2 + r14C*13CO2 + (r14C-1)*14CO2 + r14C*CO18O = 0 (14C)
        # 2r18O*CO2 + 2r18O*13CO2 + 2r18O*14CO2 + (r18O*-1)*CO18O = 0 (18O)
        # CO2 + 13CO2 + 14CO2 + CO18O = 1 (total CO2 mass balance)
        
        ## method 1: solve as matrix:
        # m =     np.array([  [r13,        -1,            0,          r13],
                            # [r14,       r14,        r14-1,          r14],
                            # [2*r18,   2*r18,        2*r18,      (r18-1)],
                            # [1,           1,            1,            1]])
        # m = m.astype('float_')
        # print(np.linalg.cond(m)) # check matrix condition number
        # n =    np.array([0, 0, 0, 1])        
        # x =    np.linalg.solve(m, n)
        # c12o_2 = x[0]
        # c13o_2 = x[1]
        # c14o_2 = x[2]
        # c12o18o = x[3]
        
        ## method 2: simultaneous equation solution:
        c13o_2 = r13*(1-r14)/(r13+1)
        c14o_2 = r14
        c12o18o = 2*r18/(1+r18)
        c12o_2 = 1 - c13o_2-c14o_2-c12o18o

        # print("Matrix vs simultaneous:") # precision check
        # print("%g\t%g" % (mc12o_2, c12o_2))
        # print("%g\t%g" % (mc13o_2, c13o_2))
        # print("%g\t%g" % (mc14o_2, c14o_2))
        # print("%g\t%g\n" % (mc12o18o, c12o18o))
               
        return c12o_2, c13o_2, c14o_2, c12o18o

    def get_gasphase_string(self):
        """Return PHREEQC input defining the object as a CO2-N2 gas phase."""
        
        # calculate 18O/O ratios in O2(g) and CO2(g)
        t_k = self.temperature + 273.15
        a_o2_h2o =  self.s.reader.get_alpha('18O','O2(g)/H2O(l)', t_k)
        a_co2_h2o = self.s.reader.get_alpha('18O','CO2(g)/H2O(l)', t_k)
        
        r18O_h2o = self.s.stnd18O * (1 + 0.001 * self.d18O)
        r18O_o2 = r18O_h2o * a_o2_h2o
        r18O_co2 = r18O_h2o * a_co2_h2o
        d18O_co2 = ( r18O_co2/self.s.stnd18O - 1 ) * 1000
        
        # get CO2 isotopologues
        iso = self.co2_isotopologues( d18O=d18O_co2 )
        
        # get O2 isotopologues
        o18o16 = 2*r18O_o2 / (1 + r18O_o2)
        o16o16 = 1 - o18o16
        
        fmt = { 'name'        :   self.name,
                'pressure'    :   PRESSURE,
                'volume'      :   self.volume,
                'temperature' :   self.temperature,
                'pp_co2'      :   PRESSURE * iso[0] * self.pCO2 * 1e-6,
                'pp_C18O'     :   PRESSURE * iso[3] * self.pCO2 * 1e-6,
                'pp_C18O2'    :   0,
                'pp_13C'      :   PRESSURE * iso[1] * self.pCO2 * 1e-6,
                'pp_13C18O'   :   0,
                'pp_13C18O2'  :   0,
                'pp_14C'      :   PRESSURE * iso[2] * self.pCO2 * 1e-6,
                'pp_14C18O'   :   0,
                'pp_14C18O2'  :   0,
                'pp_O2'       :   PRESSURE * self.O2 * o16o16,
                'pp_O18O'     :   PRESSURE * self.O2 * o18o16,
                'pp_N2'       :   PRESSURE * (1 - self.pCO2 * 1e-6 - self.O2) }        
        
        gas_string = CO2_GAS.format(**fmt)
        return gas_string
      
    def get_empty_gasphase_string(self):
        """Return an empty PHREEQC CO2-N2 gas phase string."""
        
        fmt = { 'name'        :   '',
                'pressure'    :   PRESSURE,
                'volume'      :   0,
                'temperature' :   self.s.settings['temperature'],
                'pp_co2'      :   0,
                'pp_C18O'     :   0,
                'pp_C18O2'    :   0,
                'pp_13C'      :   0,
                'pp_13C18O'   :   0,
                'pp_13C18O2'  :   0,
                'pp_14C'      :   0,
                'pp_14C18O'   :   0,
                'pp_14C18O2'  :   0,
                'pp_O2'       :   0,
                'pp_O18O'     :   0,
                'pp_N2'       :   0}

        empty_gasphase = CO2_GAS.format(**fmt)

        return empty_gasphase
   
    def mix_gases(self, Gas2, frac, **kwargs):
        """Mixes two CO2 gas phases and returns the resultant gas.
        
        Args:
            Gas2: A Gas object to be mixed with self.
            frac: The fraction (0-1) of Gas2 present in the final mixture.
            kwargs: passed to Gas constructor for new mixed Gas
        """
        
        a = self.co2_isotopologues()
        b = Gas2.co2_isotopologues()
        
        f1 = (1-frac)*self.pCO2 # total co2 from A
        f2 = frac*Gas2.pCO2     # total co2 from B
        
        mix_pco2 = f1 + f2
        mix_o2 = (1-frac)*self.O2 + frac*Gas2.O2
        m_12c = (f1*(a[0]+a[3]) + f2*(b[0]+b[3])) / mix_pco2
        m_13c = (f1*a[1] + f2*b[1]) / mix_pco2
        m_14c = (f1*a[2] + f2*b[2]) / mix_pco2
        
        mix_d13c, mix_r14c = ccu.c14_to_pmc( m_12c, m_13c, m_14c, 
                                            self.s.stnd13C, self.s.stnd14C )
        
        mix = {     'pCO2'          :   mix_pco2,
                    'd13C'          :   mix_d13c,
                    'R14C'          :   mix_r14c,
                    'O2'            :   mix_o2,
                    'temperature'   :   self.temperature,
                    'Sim'           :   self.s,
                    'name'          :   'mixed_gas'
        }
        
        kwargs.update(mix)
        return Gas(**kwargs)
   
class Simulator(object):
    """Runs Cavecalc models.
    
    A Simulator object sits at the centre of a Cavecalc model. It initalises
    IPhreeqc, executes PHREEQC input returned from Solution, Gas and 
    Carbonate objects, and handles and stores output from IPhreeqc.
    
    Simulator objects also control interaction with various utility classes
    defined in util.py.
    
    Selected Methods:
        
    """
    
    def __init__(self, settings, id=0):     
        """Initialises an IPhreeqc session.
        
        Initalising a Simulator object sets up an IPhreeqc session, parses
        the settings provided and initialises various utility objects.
        
        Args:
            settings (dict): Contains model settings. This dict should be
                generated using setter.py.
            id (int): Id number for the model being run.
        """                

        if type(settings) is not dict:
            try: 
                settings = settings.dict()
            except:
                raise TypeError("Settings arg must be of type dict.")
            
        self.settings = deepcopy(settings)         # will be updated and used
        self.settings_archive = deepcopy(settings) # will not be changed
        self.id = id
        
        # set up input log if requested
        if self.settings['phreeqc_log_file']:
            l = self.settings['phreeqc_log_file_name']
            k = l.format(self.id)
            db = self.settings['database']
            self.input_log = ccu.PhreeqcInputLog(filename = k,
                                                 dbpath = db)
        
        # set up database reader to get thermodynamic data & reactions
        self.reader = ccu.DBReader(self.settings['database'])
        
        self.count = 0
        self._bedrock_set = False
        self.string_buffer = ""
        self.desc_buffer = ""
        self._define_selected_output()
        self._parse_settings()
        
    def _define_selected_output(self):
        """Defines PHREEQC SELECTED_OUTPUT block."""
        
        fmt = { 'totals'        :   self.settings['totals'],
                'molalities'    :   self.settings['molalities'],
                'isotopes'      :   self.settings['isotopes']}
        
        if (self.settings['kinetics_mode'] == 'degas_only' or 
            self.settings['kinetics_mode'] == 'diss_only'):
            fmt['saturation_indices'] = ''
        else:
            fmt['saturation_indices'] = '\n\t-saturation_indices Calcite'
                
        self.selected_output     = SEL_OUTPUT_CORE.format(**fmt)

    def _parse_settings(self):
        """Performs any necessary settings pre-processing.
        
        1)  Calculate the initial gas (may be a mixing product)
        """
        
        ae = self.settings['atmo_exchange']
        G1 = Gas(Sim=self)
        G2 = Gas(Sim=self)
        
        G1.set_soil_gas()
        G2.set_atmosphere()
        G3 = G1.mix_gases(G2, ae)

        # Set soil water carbon & isotope budget, based on inputs
        for p in ['pCO2', 'd13C', 'R14C', 'O2']: # d13C is of CO2(g)
            if self.settings['init_'+p] == 'mix':
                self.settings['init_'+p] = getattr(G3, p)
            elif self.settings['init_'+p] == 'atm':
                self.settings['init_'+p] = getattr(G2, p)
            elif self.settings['init_'+p] == 'soil':
                self.settings['init_'+p] = getattr(G1, p)
                
    @property
    def stnd13C(self):
        if not hasattr(self, '_stnd13C'):
            self._stnd13C = self.reader.get_iso_stnd('13C')
        return self._stnd13C
    
    @property
    def stnd14C(self):
        if not hasattr(self, '_stnd14C'):
            self._stnd14C = self.reader.get_iso_stnd('14C')
        return self._stnd14C
    
    @property
    def stnd18O(self):
        if not hasattr(self, '_stnd18O'):
            self._stnd18O = self.reader.get_iso_stnd('18O')
        return self._stnd18O
    
    @property
    def stnd44Ca(self):
        if not hasattr(self, '_stnd44Ca'):
            self._stnd44Ca = self.reader.get_iso_stnd('44Ca')
        return self._stnd44Ca
       
    def _get_selected_output(self):
        """Get results from most recent Iphreeqc step. 
        
        Model parameters returned depend on SELECTED_OUTPUT. Model-essential 
        outputs are written into phreeqc_templates.py. Additional desired 
        outputs should be specified in the model input parameters (total,
        molalities, isotopes).
        
        Returns:
            A dict of outputs, as specified by SELECTED_OUTPUT.
        """
        output = self._call_iphreeqc('get_selected_output_array')
        assert len(output) > 0, "GetSelectedOutputArray returned null."
        header = output[0]
        out = {}
        for head in header:
            out[head] = []
        for row in output[1:]:
            for col, head in enumerate(header):
                out[head].append(row[col])
        
        ## correct bulk solution R14C for isotope fractionation
        # Some studies use uncorrected a14C (e.g. Genty et al. (1997)), some
        # corrected (e.g. Noronha et al. (2014)). Usually it should be
        # corrected.
        new = []
        for c13,c14 in zip(out['I_R(13C)'], out['I_R(14C)']):
            pmc = ccu.pmc_normalise(  R14C=c14,   
                                      d13C=c13, 
                                      stnd14C=self.stnd14C  )
            new.append( pmc )
        out['I_R(14C)'] = new
        return out
        
    def _output_add(self, last_only=True):
        """Add current model output to long-term storage."""
        
        if not hasattr(self, 'output'):
            self.output = {}
            
        for key in self.last_output:
            if key not in self.output.keys():
                self.output[key] = []
            if last_only:
                self.output[key].append(self.last_output[key][-1])
            else:
                for e in self.last_output[key]:
                    self.output[key].append(e)
        

    def _verify_last_output(self):
        """Check last phreeqc step ran as intended."""
        # Check phreeqc error percentage
        for value in self.last_output['pct_err']:
            if abs(value) > PHREEQC_TOLERANCE:
                raise ValueError('PHREEQC percentage error exceeds tolerance.')
    
    def _call_iphreeqc(self, command, arg=None):
        """Run a command through IPhreeqc. 
        
        Also initialises IPhreeqc if necessary. This method should not be
        called directly for reaction calls as it does not store IPhreeqc 
        output. Use ipq_buffer() and ipq_exec() instead."""
        
        ## Initialise IPhreeqc interface and database
        if not hasattr(self, 'phreeqc'):
        
            if 'phreeqpy' in MODE:
                phreeqpy_path = os.path.dirname(phreeqc_mod.__file__)
                if 'linux' in sys.platform: # linux
                    f_name = 'libiphreeqc.so.0.0.0'
                elif sys.platform == 'win32': # windows
                    f_name = 'IPhreeqc.dll'
                else: # assumed mac
                    f_name = 'libiphreeqc.0.dylib'
                try:
                    self.phreeqc = phreeqc_mod.IPhreeqc(dll_path = 
                                   os.path.join(phreeqpy_path, f_name))
                except:
                    phreeqpy_path = os.path.dirname(phreeqc_mod.__file__)             
                    e = "\nPhreeqpy failed to find the neccessary file " + \
                    "for your system (.dll/.dylib/.so). Follow instructions " + \
                    "on the phreeqpy website to compile an appropriate " + \
                    "library file and copy it to:\n%s" % os.path.join(
                    phreeqpy_path, f_name)
                    raise Exception("IPhreeqc failed to load." + e)
                    
            else: # if using API_com mode
                try:
                    self.phreeqc = Dispatch("IPhreeqcCOM.Object")
                except:
                    raise Exception("IPhreeqc failed to load. Ensure " + \
                    "IPhreeqcCOM is installed correctly.")
                    
            # IPhreeqc is now initialised. Load the database:
            self._call_iphreeqc('load_database',self.settings['database'])
         
        if 'phreeqpy' in MODE:
            c = inflection.underscore(command)
        else:
            c = inflection.camelize(command)
        f = getattr(self.phreeqc,c)
        if arg: return f(arg)
        else: return f()
    
    def ipq_buffer(self, pq_input, description=''):
        """Adds IPhreeqc input strings to the buffer.
        
        Input strings are added to the buffer and decorated with USE, SAVE and
        END statements. The buffer is used to temporarily store model input
        before it is executed by IPhreeqc. The buffer may be executed using
        the ipq_exec() method.
        
        Args:
            pq_input: A valid PHREEQC input string. Multiple strings may be
                      passed as a list of tuples [(pq_input, description), ..]. 
                      If so they will be buffered in order, with no END 
                      statements between.
            description: A short string describing the operation performed by
                         pq_input. Not necessary if pq_input is passed as a 
                         list.
        """
        
        # if multiple arguments are passed, reduce them to single strings
        if type(pq_input) is list:
            pq_statements = [a for (a,b) in pq_input]
            descriptions = [b for (a,b) in pq_input]
            pq_input = '\n'.join(pq_statements)
            description = ', '.join(descriptions)
        
        self.count += 1
        if self.count == 1:
            input = self.selected_output + pq_input + "SAVE solution 1\nEND\n"
        else:
            soln_use = "\nUSE solution %i\n" % (self.count - 1)    
            soln_save = "\nSAVE solution %i\nEND\n" % self.count
            input = soln_use + pq_input + soln_save            
            
        self.string_buffer += input
        if len(self.desc_buffer) == 0:
            self.desc_buffer = description
        else:
            self.desc_buffer += ' | ' + description

    def ipq_exec(self):
        """Run the current string buffer through IPhreeqc and log the output.
        
        Commands are added to the buffer using the ipq_buffer() method. Output 
        from the most recent call to ipq_exec() is stored in self.last_output. 
        Previous model steps are logged in self.output.
        """
        
        inp = self.string_buffer
        
        if self.settings['phreeqc_log_file']:
            self.input_log.add(inp)    #add input to log file
        self._call_iphreeqc('run_string',inp)
        
        # save outputs
        self.last_output = self._get_selected_output()
        self.last_output['step_desc'] = [self.desc_buffer]
        
        self._output_add()
        self._verify_last_output()
        
        self.desc_buffer = ""       # clear description buffer
        self.string_buffer = ""     # clear string buffer
        
    def get(self, key):
        """Look up a specified value in the current model state.
        
        Looks up the given header as a key in the self.last_output dict,
        returning the most recent value stored there. This is used to perform
        model-state-dependent calculations (e.g. isotope fractionation).
        
        Args:
            key (str): A key in self.last_output.
        Returns:
            The most recent value of this header.
        """
        data = self.last_output
        assert key in list(data.keys()), \
            ("Key %s does not exist in selected_output. Keys:\n%s" % 
            ( key, data.keys()))
        
        return data[key][-1]        # return most recent value in that entry

    def run(self):
        """Run Cavecalc!
        
        Runs a cavecalc model according to the model input parameters passed
        to the Simulator constructor. This high-level method is typically
        called immediately after initialising a Simulator object.
        
        The scriping followed depends on the kinetics_mode specified.
        
        If you want to understand the Cavecalc source code or scripting, 
        this is a good place to start.
        
        Returns:
            A dict containing step-by-step model results as specified by 
            SELECTED_OUTPUT.
        """
        
        # initialise soil water and calculate it's initial state
        water = Solution(self)
        self.ipq_buffer( *water.set_soil_chem() )
        self.ipq_exec()
        
        # initialise bedrock & calculate bedrock dissolution reaction
        bedrock = Carbonate(self)
        if self.settings['reprecip']:
            self.ipq_buffer(*water.open_diss_rxn_with_reprecip( bedrock ))
        else:
            self.ipq_buffer(*water.open_diss_rxn( bedrock ))
        self.ipq_exec()
        
        ## From here, behaviour depends on the run mode / kinetics_mode
        # kinetics mode options: 
        #   closed_system_rayleigh
        #   open_system_single_step
        #   allow_supersaturation
        #   allow_supersaturation_max
        #   degas_only
        #   dissolve_only
        
        # closed_system_rayleigh
        if self.settings['kinetics_mode'] == 'closed_system_rayleigh': # keep at SI_calcite = 0
            while water.will_it_degas(self.settings['cave_pCO2']):
                self.ipq_buffer( [water.kinetic_degas_rxn(), water.precipitate_rxn()] )
                self.ipq_exec()

        # open_system_single_step
        elif self.settings['kinetics_mode'] == 'open_system_single_step':
            self.ipq_buffer( [water.eq_degas_rxn(), water.precipitate_rxn()] )
            self.ipq_exec()
                
        # allow_supersaturation
        elif self.settings['kinetics_mode'] == 'allow_supersaturation':
            while water.will_it_degas(self.settings['cave_pCO2']):
                while self.get('si_Calcite') < self.settings['calcite_sat_limit']:
                    self.ipq_buffer( *water.kinetic_degas_rxn() )
                    self.ipq_buffer( WATER_EQUILIBRATE, 'DIC equilibrate' )
                    self.ipq_exec()
                self.ipq_buffer( *water.precipitate_rxn() )
                self.ipq_exec()

        # allow_supersaturation_max
        elif self.settings['kinetics_mode'] == 'allow_supersaturation_max':
            while water.will_it_degas(self.settings['cave_pCO2']):
                while water.will_it_degas(self.settings['cave_pCO2']):
                    self.ipq_buffer( *water.kinetic_degas_rxn() )
                    self.ipq_buffer( WATER_EQUILIBRATE, 'DIC equilibrate' )
                    self.ipq_exec()
                self.ipq_buffer( *water.precipitate_rxn() )
                self.ipq_exec()
                
        # degas_only
        elif self.settings['kinetics_mode'] == 'degas_only':
            while water.will_it_degas(self.settings['cave_pCO2']):
                self.ipq_buffer( *water.kinetic_degas_rxn() )
                self.ipq_buffer( WATER_EQUILIBRATE, 'DIC equilibrate' )
                self.ipq_exec()
                
       # diss_only
        elif self.settings['kinetics_mode'] == 'diss_only': # dissolve bedrock only
                pass
                
        ccu.PostProcessor(self) # perform offline calculations and processing
        
        return self.output
        
    def save_results(self, format='.pkl', filename='output'):
        """Save model output to specified file and format.
        
        Files will be saved in the current directory.
        
        This method is provided for low-level access. Generally it is advisable
        to run Cavecalc via forward_models.py and use the more powerful saving
        functions of that module.
        
        Args:
            format (str): '.csv', '.pkl' or '.mat'
            directory (str): Save directory
            filename (str): Save file name without extension. Not used for pkl.
        """
                
        if 'csv' in format:
            ccu.save_csv(self.output, '%s.csv' % filename)
        elif 'pkl' in format:
            ccu.save_pkl([self.output], 'results.pkl')
            ccu.save_pkl([self.settings], 'settings.pkl')
        elif 'mat' in format:
            ccu.save_mat(self.output, '%s.mat' % filename)
        else:
            raise IOError("Unrecognised output format.")