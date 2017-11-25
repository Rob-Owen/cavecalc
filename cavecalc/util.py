"""Cavecalc module for utility functions and classes.

Contains code to perform file IO, model output post-processing and various
other helpful bits.
"""

import os
import datetime
import math
import pickle
import csv
import re
from collections import OrderedDict
from copy import copy
import numpy as np
import scipy.io as sio
import cavecalc.data

class DBReader(object):
    """Reads data from the PHREEQC database file.
    
    DBReader objects are used to extract out database information (e.g. 
    partition coefficients and fractionation factors) that may be needed in
    geochemical calculations coded in caves.py. Thermodynamic parameters are
    not hard-coded into Cavecalc python files.
    
    DBReader objects cache use regular expressions to parse the database file 
    and cache data extracted to improve performance on repeated requests.
    
    Selected Methods:
        get_k_values - get phase definition thermodynamic data
        get_1000lnalpha - get an isotope fractionation factor
        get_alpha - get an isotope fractionation factor
        get_iso_stnd - get the absolute isotope ratio in a standard
    """
    
    def __init__(self, database):
        """Initialise the object.
        
        Args:
            database (str): location of phreeqc database filename.
        """
        
        db_dir = os.path.dirname(cavecalc.data.__file__)
        self.db = os.path.join(db_dir,database)
        
        self.alphas = {}
        self.thermo = {}
        
    def _cc(self, line):
        return line.split('#')[0]
        
    def _database_eval(self, equation, temperature): 
        """Evaluate a database expression for a fractionation factor or
        partition coefficent.
        
        Args:
            equation (list): A list of coefficents for a PHREEQC thermodynamic
                            parameter.
            temperature: Temperature (degrees K) at which to evaluate the
                            equation.
        Returns:
            Value of the expression at the given temperature
        """
        T = temperature
        a = equation
        
        coeffs = [1, T, 1/T, math.log10(T), 1/(T*T), T*T]        
        evaluation = 0
        for i in range(len(a)):
            evaluation += float(a[i])*float(coeffs[i])   # compute 1000ln_alpha
        return evaluation
            
    def _phase_lookup(self, reactants_list):
        
        r = r' \+ '.join(reactants_list)
        rc1 = re.compile(r"\s*?{}\s*=".format(r))
            
        out = []
        
        with iter(open(self.db, 'r')) as f:
            for line in f:
                if re.match(rc1, line):
                    line = next(f)
                    while line.strip() not in ('', '#'):
                        a = self._cc(line.strip()).split()
                        out.append(copy(a))
                        line = next(f)
                    return [o for o in out if o]
        raise Exception("No PHASES entry matched: %s" % r)
        
    def _ne_lookup(self, isotope, species):
        
        def escape_brackets(string):
            s = string.replace(r'(', r'\(')
            s = s.replace(r')',r'\)')
            s = s.replace(r'\\)',r'\)')
            s = s.replace(r'\\(',r'\(')
            s = s.replace(r'\[',r'[')
            s = s.replace(r'\]',r']')
            return s
        
        isotope = escape_brackets(isotope)
        species = escape_brackets(species)
            
        r1 = r"\s*Log_alpha_{}_{}\s*".format(isotope, species)
        rc1 = re.compile(r1)
        out = []
        
        with iter(open(self.db, 'r')) as f:
            for line in f:
                if re.match(rc1, line):
                    line = next(f)
                    while line.strip() not in ('', '#'):
                        a = self._cc(line.strip()).split()
                        out.append(copy(a))
                        line = next(f)
                    return [o for o in out if o]
        raise Exception("No NAMED EXPRESSIONS entry matched: %s" % r1)
                
    def get_k_values(self, reactants_list, temperature=298.15):
        """Get thermodynamic data for a specified phase.
        
        Looks up the reaction specified by 'reactants_list' in the database
        and returns all thermodynamic data found.
        
        Args:
            reactants_list: A list of strings. Each string defines a species
                such that the species listed, all reacted together, specify a
                dissolution reaction (i.e. PHASES defintiion) in the database
                file.
            temperature: The temperature to evaluate any analytical expressions
                at.
                
        Returns:
            A dict containing data found:
            log_k: log_k for the reaction, as defined in the database.
            delta_h: delta_h for the reaction, as defined in the database.
            gamma: gamma for the reaction, as defined in the database.
            analytic_line: string containing analytical K expression constants
            analytic: The evaluation of analytic_line at the specified
                temperature.
                
            Entries not present in the database are returned as None.
        """
        
        if temperature < 50:
            print("Warning! Temperature (K) given as %i. Converting to %i" %
                    (temperature, temperature+273.15))
            temperature = temperature + 273.15
    
        # If possible return cached values
        chk = reactants_list.copy()
        chk.append(str(temperature))
        chk = ' '.join(chk)
        if chk in self.thermo:
            return self.thermo[chk]
            
        thermo = {}
        data = self._phase_lookup(reactants_list)
        for a in data:
            if a[0] == 'log_k':
                thermo[a[0]] = a[1]
            elif a[0] == 'delta_h':
                thermo[a[0]] = a[1]
            elif a[0] == '-gamma':
                thermo[a[0]] = a[1]
            elif a[0] == '-analytic':
                thermo['analytic_line'] = ' '.join(a)
                thermo['analytic_value'] = self._database_eval( a[1:], 
                                                                temperature )
        
        self.thermo[chk] = thermo
        return thermo
                
    def get_1000lnalpha(self, isotope, species, temperature=298.15):
        """Get an isotope fractionation factor from the database.
        
        Returns 1000lnalpha for the specified reaction.
        
        Args:
            isotope (str): the isotope fractionating (e.g. '13C')
            species (str): the species pair involved in the fractionation
                (eg. r'CO2(g)/CO2(aq)')
            temperature: The temperature (degrees K) to calculate 1000lnalpha
                at.
        Returns:
            1000lnalpha for the reaction.
        """
        
        if temperature < 50:
            print("Warning! Temperature (K) given as %i. Converting to %i" %
                    (temperature, temperature+273.15))
            temperature = temperature + 273.15
            
        # Check for cached value
        chk = isotope + species + str(temperature)
        if chk in self.alphas:
            return self.alphas[chk]
        
        
        value = 0
        tmp = self._ne_lookup(isotope, species)
        
        for a in tmp:
            if a[0] == '-add_logk':
                b = a[1].split('_')
                value += self.get_1000lnalpha(b[2], b[3], temperature)*int(a[2])
            elif a[0] == '-ln_alpha1000':
                value += self._database_eval(a[1:], temperature)                
        self.alphas[chk] = value   
        return value
   
    def get_alpha(self, isotope, species, temperature=298.15):
        """Get an isotope fractionation factor from the database.
        
        Returns alpha for the specified reaction. Arguments are the same as
        get_1000lnalpha.
        
        Args:
            isotope (str): the isotope fractionating (e.g. '13C')
            species (str): the species pair involved in the fractionation
                (eg. r'CO2(g)/CO2(aq)')
            temperature: The temperature (degrees K) to calculate 1000lnalpha
                at.
        Returns:
            alpha for the reaction.
        """
    
        ln1000a = self.get_1000lnalpha(isotope, species, temperature=298.15)
        return math.exp(0.001*ln1000a)
    
    def get_iso_stnd(self, isotope):
        """Get the absolute isotope ratio in the standard.
        
        Looks up the absolute mole fraction of the specified isotope present
        in the standard. The standard used, and it's composition, are defined
        in the phreeqc database.
        
        Args:
            isotope (str): The isotope desired (e.g. '13C')
        Returns:
            The mole fraction of 'isotope' present in the standard.
        """
        
        a = None
        pat = "([ \t]*?)-isotope([ \t]*?)\[{!s}\][ \t]"
        c = re.compile(pat.format(isotope))
        with open(self.db, 'r') as f:
            for line in f:
                if re.match(c, line):
                    a = line
                    
        if a is None:
            raise Exception("No isotope standard found.")
            
        b = a.split()
        return float(b[3])
            
class PhreeqcInputLog(object):
    """Logs IPhreeqc input to a text file.

    PhreeqcInputLog is a utility class used by caves.Simulator objects to log
    iphreeqc input strings to a .phr text file. The resulting file is useful 
    for understanding how the code runs and debugging failed models. The log 
    file is also valid phreeqc input and can be run directly as a phreeqc 
    script.
    """
    
    def __init__(self, filename, dbpath):
        """Initialise the object and create a log file.
        
        Creates a .phr log file at the specified location. The first few lines
        are initialised with the date, time and full path to the database used.
        
        Args:
            filename (str): Location to write the log file.
            dbpath (str): Path to the database file. Perferably full path.
        """
        
        self.filename = filename
        self.pq_input = open(self.filename,'w')
        self.pq_input.truncate()
        self._preamble()
        
        self.pq_input.write("\nDATABASE %s" % dbpath)
        self.pq_input.close()
        
    def _preamble(self):
        """Write log file preamble."""
        
        now = datetime.datetime.now()
        line1 = "#\tIPhreeqc input log.\n"
        line2 = "#\tDate:\t%i-%i-%i\n" % \
            (now.year, now.month, now.day)
        line3 = "#\tTime:\t%i:%i:%i\n" % \
            (now.hour, now.minute, now.second)
            
        self.pq_input.write(line1 + '\n')
        self.pq_input.write(line2)
        self.pq_input.write(line3)
        
    def _buffer(self):
        """Write break into log file for readability."""
        
        self.pq_input.write("\n\n" + '#' + '-'*20 + "\n\n")
    
    def add(self, string):
        """Write a string to log file.
        
        Args:
            string (str): PHREEQC input text to be added to log file.
        """
        self.pq_input = open(self.filename,'a')
        self._buffer()
        self.pq_input.write(string)
        self.pq_input.close()

class PostProcessor(object):
    """Performs offline calculations and formatting of model results.
    
    Post-processes model output from a caves.Simulator object, adding 
    parameters calculated offline and making the output more readable.
    """
    
    def __init__(self, Simulator):
        """Run the post-processor.
        
        Does the following:
            - Adds f_ca (fraction of Ca remaining) and f_c (fraction of c 
              remaining) parameters to the model output.
            - Calculates X/Ca ratios from mole outputs (Mg, Sr, Ba)
            - Rename some outputs to be more readable
            - Remove nonsense radiocarbon calculations from the degassing loop
              results.
            - Formats some PHREEQC-output parameter names to make them more
              readable.
        
        Args:
            Simulator: The Simulator object to operate on.
        """
        
        self.s = Simulator
        self.calculate_f()
        self.calculate_XCa()
        self.tidy()
        self.calculate_radiocarbon()
        self.set_none()

    def calculate_f(self):
        """Calculate the fractions of carbon and calcium remaining at each
        model step.
        
        Calculates f_c (fraction of C) and f_ca (fraction of Ca) remaining in
        the solution at each model step relative to the amount present
        immediately following bedrock dissolution. Calculated values are added 
        to the Simulator output dict.
        """

        init = output_filter(self.s.output, 'step_desc', 'dissolve')
        init_ca = init['Ca(mol/kgw)'][0]
        init_c = init['C(mol/kgw)'][0]
        ca = self.s.output['Ca(mol/kgw)']
        c = self.s.output['C(mol/kgw)']
        self.s.output['f_ca'] = [x / init_ca for x in ca]
        self.s.output['f_c'] = [x / init_c for x in c]
        
    def calculate_XCa(self):
        """Calculate X/Ca ratios of solution and precipitates.
        
        X/Ca is calculated at each time step using a Rayleigh distillation
        model. X/Ca ratios for solution and any precipitate are added to the
        Simulator output dict.
        
        X includes Mg, Sr and Ba.
        """
        
        a = self.s.output
        
        trace_elements = ['Ba', 'Sr', 'Mg']
        dissolved_ratios = {k : [] for k in trace_elements}
        precipitate_ratios = {k : [] for k in trace_elements}
        
        for i, desc in enumerate(a['step_desc']):
            w = a['mass_H2O']
            ca = a['Ca(mol/kgw)']
            dissolved_ca = ca[i] * w[i]
            
            if 'CaCO3_precipitation' in desc:
                solid_ca = dissolved_ca - ca[i-1] * w[i-1]
                
            for x in trace_elements:
                x_dissolved = a[x+'(mol/kgw)'][i] * w[i]
                try: 
                    dissolved_ratios[x].append( x_dissolved / dissolved_ca )
                except ZeroDivisionError:
                    dissolved_ratios[x].append( 0 )
                
                # precipitated trace element = amount lost from solution
                if 'CaCO3_precipitation' in desc:
                    x_precip = x_dissolved - a[x+'(mol/kgw)'][i-1] * w[i-1]
                    try: 
                        precipitate_ratios[x].append( x_precip / solid_ca )
                    except ZeroDivisionError:
                        precipitate_ratios[x].append( 0 )
                else:
                        precipitate_ratios[x].append( 0 )
         
        for x in trace_elements:
            self.s.output[x+'/Ca(mol/mol)'] = dissolved_ratios[x]
            self.s.output[x+'/Ca(mol/mol)_Calcite'] = precipitate_ratios[x]
                
    def calculate_radiocarbon(self):
        """Sets all post-dissolution R14C values to be constant.
        
        Radiocarbon distribution is fully calculated during the bedrock 
        dissolution step, but is not included in subsequent degassing /
        precipitation steps for three reasons:
            - It increases computation time.
            - When model steps are very small, 14C may cause non-convergence.
            - d13C-corrected R14C is not expected to change during these 
                processes.
            
        To avoid confusion, spurious post-dissolution R14C values reported by
        IPhreeqc are set equal to the post-dissolution value.
        
        This method also adds a new column ('DCP') to the output. Note that
        the DCP calculation assumes modern = 100 pMC.
        """

        for key in self.s.output.keys():
            if 'R14C' in key:
                pmc = None
                for (i, value) in enumerate(self.s.output[key]):
                    if 'bedrock' in self.s.output['step_desc'][i]:
                        pmc = self.s.output[key][i]
                    elif pmc is not None:
                        self.s.output[key][i] = pmc
                        
        atmo_a14c_init = 100 # placeholder value (pMC)
        self.s.output['DCP'] = [ (1-v/atmo_a14c_init)*100 for v in 
                                 self.s.output['R14C'] ]
                        
    def tidy(self):
        """Rename PHREEQC variable names with more reader-friendly syntax. 
        
        Also removes a couple of unwanted outputs.
        """
        o = self.s.output
        rep_keys = {}
        
        # remove unwanted outputs (not of interest to most users)
        o.pop('soln')
        o.pop('mass_H2O')
        o.pop('pct_err')
        o.pop('temp(C)')
        o.pop('I_R(14C)_CO2(aq)')
        
        # rename some outputs for readability
        for key, data in o.items():
            if key[0:3] == 'I_R': # rename isotopes
                iso = key[key.find('(')+1:key.find(')')]
                if iso == '14C':
                    new_key = 'R' + iso + key[key.find(')')+1:]
                else:
                    new_key = 'd' + iso + key[key.find(')')+1:]
                rep_keys[key] = new_key
            elif key[0:2] == 'm_': # rename molalities
                new_key = key[2:]
                rep_keys[key] = new_key
            elif key[0:2] == 's_':
                new_key = 'moles_' + key[2:]
                rep_keys[key] = new_key
        
        for k1, k2 in rep_keys.items():
            o[k2] = o.pop(k1)
            
    def set_none(self):
        """Set default paramter returns (e.g. d18O = -999) to None."""
        
        for k, lst in self.s.output.items():
            new_lst = []
            for e in lst:
                if not(isinstance(e, str)) and e <= -999:
                    new_lst.append(None)
                else:
                    new_lst.append(e)
            self.s.output[k] = new_lst
       
# Begin radiocarbon calculation functions
def pmc(C14, d13C, stnd14C=1.175887709e-12):
    """Converts 14C/C absolute ratio to a d13C-corrected pMC value.
    
    Args:
        C14 : 14C/C ratio (absolute)
        d13C : V-PDB d13C values
        stnd14C : The 14C/C ratio in the standard
    Returns:
        the standard-normalised radiocarbon value in pMC.
    """
    
    return C14 / stnd14C * 100 * pow(0.975/(1+0.001*d13C),2)
       
def pmc_2_c14(R14C, d13C, stnd14C=1.175887709e-12):
    """Converts a d13C-corrected pMC value to a 14C/C absolute ratio.
    
    Conversion follows Stuvier & Pollach (1977).
    
    Args:
        R14C: A radiocarbon value in d13C corrected percent modern carbon
              (pmc).
        d13C: The corresponding d13C value.
        stnd14C: The standard 14C/C molar ratio. Default is 
                 1.175887709e-12.
    Returns:
        14C/C ratio
    """
    
    return stnd14C * 0.01 * R14C * pow((1+0.001*d13C)/0.975,2)
       
def c14_to_pmc(C12, C13, C14, stnd13C=0.0111802, stnd14C=1.175887709e-12):
    """Returns isotope ratios given relative c isotope abundance data.
    
    Args:
        C12: Mole fraction 12C
        C13: Mole fraction 13C
        C14: Mole fraction 14C
        stnd13C: 13C/12C ratio in the standard. Default is 0.0111802 (VPDB)
        stnd13C: 14C/C ratio in the standard. Default is 1.175887709e-12
        
    Returns:
        d13C, R14C (pMC, d13C-corrected)
    """
    
    d13C = ((C13/C12)/stnd13C - 1) * 1000
    R14C = pmc(C14, d13C, stnd14C)
    return d13C, R14C
        
def pmc_normalise(R14C, d13C, stnd14C=1.175887709e-12):
    """Convert a 14C/C pmc-normalised ratio (returned by PHREEQC) input a 
    'proper' d13C corrected pMC ratiocarbon ratio."""
    
    true_ratio = 0.01 * R14C * stnd14C
    return pmc(true_ratio, d13C, stnd14C)
    
def pmc_denormalise(pMC, d13C, stnd14C=1.175887709e-12):
    """Convert a d13C-corrected R14C (pMC) to a normalised 14C/C ratio, as used
    by PHREEQC."""
    
    c14_ratio = pmc_2_c14(pMC, d13C, stnd14C)
    return 100 * c14_ratio / stnd14C
    
# begin helper function definitions    
def output_filter(src, key, value):
    """Filter lists contained in a dict by the value in one of the lists.
    
    The function takes a dict of lists, where all lists are of equal length. It
    returns a new list of dicts where the lists are shortened: only list 
    entries at indices that meet specified criteria are copied.
    
    Args:
        src: Dictionary to filter. Each entry should be a list of equal length.
        key: The key in src to filter the lists by
        value: value in src[key] to include in output.
    Returns:
        A filtered version of src.
    """
    
    if isinstance(value, str):
        f = lambda v : value in v
    else:
        f = lambda v : value == v
        
    inds = [i for i,a in enumerate(src[key]) if f(a)]
    
    o = {}
    for k, v in src.items():
        o[k] = [a for i,a in enumerate(v) if i in inds]
    return o
    
def matlab_header_parse(dictionary):
    """Remove illegal characters from dictionary keys.
    
    To transfer data to matlab, certain characters are not allowed in field / 
    array names (e.g. brackets, hypens). This function removes these characters
    in preparation for writing dicts (of numpy arrays) to a .mat file.
    
    Args:
        dictionary: A dict
    Returns:
        A dict with modified key names
    """
    
    b = {}
    for k in dictionary:
        new_k = k.replace(')', '') # because matlab is fussy
        new_k = new_k.replace('(','_') # about variable names
        new_k = new_k.replace('-','') # and does not allow these
        new_k = new_k.replace('/','') # characters
        new_k = new_k.replace('[','')
        new_k = new_k.replace(']','')    
        
        b[new_k] = dictionary[k]
    return b
    
def numpify(dictionary):
    """Prepare a dict of lists for writing to a .mat file.
    
    Convert the dict of lists to a dict of numpy arrays. Dict keys are edited
    if they contain matlab-illegal characters. Lists in the dict are converted
    to numpy arrays.
    
    Args:
        dictionary: A dict of lists. Each list should be composed of a single
            type.
    Returns:
        A modified dict, ready for writing to a .mat file.    
    """
    
    a = matlab_header_parse(dictionary)
    b = {}
    for k in a:
        if k == "step_desc":
            b[k] = np.asarray( a[k], order='F' )
        else:
            b[k] = np.asarray( a[k], order='F' )    
    return b
            
def save_mat(dict_of_lists, filename):
    """Save data to a .mat file for use with Matlab/Octave.
    
    Takes a dict of lists (e.g. model results) and saves them to a .mat file.
    Data are prepared for saving using the numpify() function.
    
    Args:
        dict_of_lists: Data to be saved for Matlab use.
    """
        
    out = numpify(dict_of_lists)
    sio.savemat( filename, out )

def save_pkl(data, filename):
    """Save data to a .pkl file for use with Python.
    
    Args:
        data: a pickle-able Python object.
        filename: Output file name/location.
    """
    
    with open(filename,'wb') as f:
        pickle.dump( data, f)
            
def save_csv(dictionary, filename):
    """Save data to a .csv file.

    For saving data to view in for use with other programs (e.g. Excel). Data
    should be provided as a dict of lists. Dict keys give column headers and
    lists give column values. All lists must be of equal length.
    
    In the resulting .csv, columnns are arranged alphabetically by header.
    
    Args:
        dictionary: The dictionary to write to file.
    """
    
    r = OrderedDict(sorted(dictionary.items()))
    a, b = zip(*[(k,v) for (k,v) in r.items()])
    c = zip(*b)
    
    with open(filename,'w', newline='') as f:
        writer = csv.writer(f,dialect='excel')
        writer.writerow(a)
        for row in c:
            writer.writerow(row)