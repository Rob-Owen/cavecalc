from cavecalc.configuration import RunConfig
from cavecalc.util import output_filter
from typing import Dict, List, Iterable
from copy import deepcopy

def calculate_f(output: Dict[str, List]) -> None:
    """Calculate the fractions of carbon and calcium remaining at each
    model step.
    
    Calculates f_c (fraction of C) and f_ca (fraction of Ca) remaining in
    the solution at each model step relative to the amount present
    immediately following bedrock dissolution. Calculated values are added 
    to the Simulator output dict.
    """

    init = self.filter('step_desc', 'dissolve')
    init_ca = init['Ca(mol/kgw)'][0]
    init_c = init['C(mol/kgw)'][0]
    output['f_ca'] = [x / init_ca for x in output['Ca(mol/kgw)']]
    output['f_c'] = [x / init_c for x in output['C(mol/kgw)']]
    
def calculate_XCa(output: Dict[str, List]) -> None:
    """Calculate X/Ca ratios of solution and precipitates.
    
    X/Ca is calculated at each time step using a Rayleigh distillation
    model. X/Ca ratios for solution and any precipitate are added to the
    Simulator output dict.
    
    X includes Mg, Sr and Ba.
    """
    
    a = output
    
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
        output[x+'/Ca(mol/mol)'] = dissolved_ratios[x]
        output[x+'/Ca(mol/mol)_Calcite'] = precipitate_ratios[x]
            
def calculate_radiocarbon(output: Dict[str, List]) -> None:
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

    for key in output.keys():
        if 'R14C' in key:
            pmc = None
            for (i, value) in enumerate(output[key]):
                if 'bedrock' in output['step_desc'][i]:
                    pmc = output[key][i]
                elif pmc is not None:
                    output[key][i] = pmc
                    
    atmo_a14c_init = 100 # placeholder value (pMC)
    output['DCP'] = [ (1-v/atmo_a14c_init)*100 for v in 
                                output['R14C'] ]
                    
def tidy(output: Dict[str, List]) -> None:
    """Rename PHREEQC variable names with more reader-friendly syntax. 
    
    Also removes a couple of unwanted outputs.
    """
    o = output
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

def set_none(output: Dict[str, List]) -> None:
    """Set default paramter returns (e.g. d18O = -999) to None."""
    def find_defaults(vals: Iterable) -> List[int]:
        for i, v in enumerate(vals):
            if not(isinstance(e, str)) and e <= -999:
                yield i
    
    for lst in output.values():
        lst[find_defaults(lst)] = None

class Results:
    self.run_config: RunConfig

    def add_output_step(self, step) -> None:
        if not hasattr(self, 'output'):
            self.results = {}

        for key in self.last_output:
            if key not in self.results.keys():
                self.results[key] = []
            if last_only:
                self.results[key].append(self.last_output[key][-1])
            else:
                for e in self.last_output[key]:
                    self.results[key].append(e)

    
    def get_current_value(self, key: str):
        return self.results[key][-1]

    def finalise() -> Dict[str, List]:
        # perform post-processing
        output = deepcopy(self.results)

        calculate_f(output)
        calculate_XCa(output)
        tidy(output)
        calculate_radiocarbon(output)
        set_none(output)
        return output


    def _filter(self, key, value) -> Dict[List]:
        """Returns a filtered copy of output.
        """

        if isinstance(value, str):
            f = lambda v : value in v
        else:
            f = lambda v : value == v
            
        inds = [i for i,a in enumerate(self.results([key]) if f(a)]
        o = {}
        for k, v in self.results.items():
            o[k] = [a for i,a in enumerate(v) if i in inds]
        return o

