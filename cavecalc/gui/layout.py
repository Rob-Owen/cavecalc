"""Codes locations and type of each parameter in the GUI input panel.

    0: Hidden
    10: Geochem: Second Gas End-Member
    11: Geochem: Soil End-Member
    12: Geochem: Soil gas mixing
    13: Geochem: Bedrock chemistry
    14: Geochem: Bedrock dissolution conditions
    15: Geochem: General (e.g. temp)
    16: Geochem: Soil Metals
    17: Geochem: Cave Air
    2: Scripting options
    3: Additional PHREEQC Output
    4: IO Settings
    
Type:
    A: Text-with-range
    B: Text without range
    C: Options Menu
    D: Check button
    E: File browser - load
    F: File browser - select directory
    """

database =        (4, 'E')
phreeqc_log_file =     (4, 'D')
phreeqc_log_file_name = (0, 'B') # no reason to change
temperature =          (15, 'A')
co2_decrement =        (2, 'A')
calcite_sat_limit =    (2, 'A')
bedrock =              (14, 'A')
bedrock_mineral =      (13, 'A')
bedrock_pyrite =       (14, 'A')
bedrock_d44Ca =        (13, 'A')
bedrock_d13C =         (13, 'A')
bedrock_d18O =         (13, 'A')
bedrock_MgCa =         (13, 'A')
bedrock_SrCa =         (13, 'A')
bedrock_BaCa =         (13, 'A')
atmo_exchange =		(12, 'A')
gas_volume =		(14, 'A')
atm_pCO2 =		    (10, 'A')
atm_d13C =		(10, 'A')
atm_R14C =		(10, 'A')
atm_d18O =		(15, 'A')
init_pCO2 =		(12, 'A')
init_d13C =		(12, 'A')
init_R14C =		(12, 'A')
init_solution_d13c = (0, 'A')    # not exposed in GUI useage
soil_pH =       (0, 'A')        # pH is charge-balanced - initial estimate is unimportant
soil_pCO2 =     (11, 'A')
soil_Ca =       (16, 'A')
soil_Mg =       (16, 'A')
soil_Sr =       (16, 'A')
soil_Ba =       (16, 'A')
soil_d13C =     (11, 'A')
soil_R14C =     (11, 'A')
soil_d44Ca =    (16, 'A')
cave_O2 =       (17, 'A')
cave_pCO2 =		(17, 'A')
cave_d13C =		(17, 'A')
cave_R14C =		(17, 'A')
cave_d18O =		(17, 'A')
cave_air_volume = (17, 'A')
kinetics_mode = (2, 'C')
reprecip =      (14, 'D')
totals =        (3, 'B')
molalities =    (3, 'B')
isotopes =      (3, 'B')
out_dir =       (4, 'F')
soil_O2 =       (11, 'A')
atm_O2 =        (10, 'A')
init_O2 =       (12, 'A')
