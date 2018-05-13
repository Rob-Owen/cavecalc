"""Default Cavecalc Model Input Parameters"""


# Atmospheric gas end-member
atm_O2   =      0.21 # Atmospheric O2 (%)
atm_d18O =		-10	 # Rainfall d18O (per mil VSMOW)
atm_pCO2 =		270	 # Atmospheric pCO2 (ppmv)
atm_d13C =		-7   # Atmospheric CO2 d13C (per mil VPDB)
atm_R14C =		100  # Atmospheric CO2 R14C (pmc)

# Soil gas end-member
soil_O2   =            0        # end-member soil gas O2 (%)
soil_R14C =            100      # end-member soil gas R14C (pmc)
soil_d13C =            -25	    # end-member soil gas d13C (per mil VPDB)
soil_pCO2 =            20000    # end-member soil gas pCO2 (ppmv)

# Soil Gas Mixing
atmo_exchange =		0	# Atmospheric contribution to soil gas CO2 (0-1)
init_O2   =     'mix'   # actual soil gas O2 content (%)
init_R14C =		'mix'	# actual soil gas R14C* (pmc)
init_d13C =		'mix'	# actual soil gas d13C* (per mil VPDB)
init_pCO2 =		'mix'	# actual soil gas pCO2* (ppmv). 'atm' to equilibrate with atmospheric, 'soil' to equilibrate with soil_pCO2 (after atmospheric mixing)
init_solution_d13c = 'n/a'

# Soil Metals
soil_Ba =              0    	# soil water Ba concentration (mmol/kg water)
soil_Ca =              0    	# soil water Ca concentration (mmol/kg water)
soil_Mg =              0    	# soil water Mg concentration (mmol/kg water)
soil_Sr =              0    	# soil water Sr concentration (mmol/kg water)
soil_pH =              4        # initial soil pH estimate - not important
soil_d44Ca =           0        # soil water Ca d44Ca (per mil 915a)

# Bedrock Chemistry
bedrock_BaCa =         0            # bedrock Ba/Ca (mmol/mol)
bedrock_MgCa =         0            # bedrock Mg/Ca (mmol/mol)
bedrock_SrCa =         0            # bedrock Sr/Ca (mmol/mol)
bedrock_d13C =         0            # bedrock d13C (per mil VPDB)
bedrock_d18O =         0            # bedrock d18O (per mil VSMOW)
bedrock_d44Ca =        0            # Bedrock d44Ca (per mil SRM 915a)
bedrock_mineral =      'Calcite'    # Bedrock mineralogy ('Calcite' or 'Dolomite')

# Bedrock Dissolution Conditions
bedrock =              10   # moles bedrock (defaults to excess)
bedrock_pyrite =       0	        # moles of pyrite present during bedrock dissolution
gas_volume =		   10	# Volume of soil gas present during bedrock dissolution (L/kg water)
reprecip =             False    # Allow calcite re-precipitation during bedrock dissolution

# Cave Air
cave_O2 =              0.21    
cave_pCO2 =		       1000
cave_d13C =		       -10
cave_R14C =		       100
cave_d18O =		       0
cave_air_volume =      0

# General
temperature =          20     # Temperature in degrees celsius
kinetics_mode =        'multi_step_degassing'    # Specifies how to run the model (see types_and_limits.py for options)

# Scripting Options
co2_decrement =        0.5    # Fraction of CO2(aq) removed on each degassing step
calcite_sat_limit =    1      # Only used if kinetics_mode = 'ss'. CaCO3 only precipitates when saturation index exceeds this value. 

# Additional PHREEQC output
# SELECTED_OUTPUT (model-essential outputs are set in phreeqc_templates.py)
isotopes =   """R(44Ca) R(18O) R(13C) R(18O)_HCO3- R(13C)_HCO3- R(18O)_CO3-2
                R(13C)_CO3-2 R(44Ca)_Calcite R(18O)_Calcite R(13C)_Calcite"""
molalities =   "HCO3-	CO3-2"
totals =       ""

# File IO Settings
phreeqc_log_file =      False           # If selected, a .phr file will be included in the output.
out_dir =               ""              # Save location
database =              'oxotope.dat'   # database file
phreeqc_log_file_name = 'log_{}.phr'    # The .phr file name (if saved). {} will be replaced by the model number. 




               
