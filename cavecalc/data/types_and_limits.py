"""Encodes accepted values and ranges for model input parameters.

This data is used by the setter.SettingsObject.validate_entry() method. To
verify input parameter types and values.

     'str'      :   x must be a string
     (a, b)     :   x must be numeric, with x > a and x < b. a or b == None 
                    means there is no lower or upper bound, respectively.
     (a, b, ...):   As above. If x is not numeric it must match an entry in ...
     True       :   x must be boolean
     [...]      :   x must match a entry in the list    
"""

database =        'str'
phreeqc_log_file =     True
phreeqc_log_file_name = 'str'
out_dir =              'str'
temperature =          (-273.15, None)
co2_decrement =        (0, 1)
calcite_sat_limit =    (0, None)
bedrock =              (0, None)
bedrock_mineral =      ['Calcite', 'Dolomite']
bedrock_pyrite =       (0, None)
bedrock_d44Ca =        (None, None)
bedrock_d13C =         (None, None)
bedrock_d18O =         (None, None)
bedrock_MgCa =         (0, None)
bedrock_SrCa =         (0, None)
bedrock_BaCa =         (0, None)
atmo_exchange =		(0, 1)
gas_volume =		(0, None)
atm_O2 =        (0, 100)    
atm_pCO2 =		(0, None)
atm_d13C =		(None, None)
atm_R14C =		(0, None)
atm_d18O =		(None, None)
cave_O2 =        (0, 100)    
cave_pCO2 =		(0, None)
cave_d13C =		(None, None)
cave_R14C =		(0, None)
cave_d18O =		(None, None)
cave_air_volume = (0, None)
init_O2 =       (0, 100, 'soil', 'atm', 'mix')
init_pCO2 =		(0, None, 'soil', 'atm', 'mix')
init_d13C =		(None, None, 'soil', 'atm', 'mix')
init_R14C =		(0, None, 'soil', 'atm', 'mix')
init_solution_d13c = (None, None, 'n/a')
soil_O2 =       (0, 100)
soil_pH =       (0, 14)
soil_pCO2 =     (0, None)
soil_Ca =       (0, None)
soil_Mg =       (0, None)
soil_Sr =       (0, None)
soil_Ba =       (0, None)
soil_d13C =     (None, None)
soil_R14C =     (0, None)
soil_d44Ca =    (None, None)
kinetics_mode = ['closed_system_rayleigh', 'open_system_single_step', 'allow_supersaturation', 'allow_supersaturation_max', 'degas_only', 'diss_only']
reprecip =      True
totals =        'str'
molalities =    'str'
isotopes =      'str'