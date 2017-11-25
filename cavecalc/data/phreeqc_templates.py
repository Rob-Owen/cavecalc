"""PHREEQC format strings - used by caves.py"""

SEL_OUTPUT_CORE = """
SELECTED_OUTPUT
    -high_precision     true
    -reset              false
    -simulation         false
    -state              false
    -solution           true
    -distance           false
    -time               false
    -step               false
    -ph                 true
    -pe                 false
    -reaction           false
    -temperature        true
    -alkalinity         false
    -ionic_strength     false
    -water              true
    -charge_balance     false
    -percent_error      true{saturation_indices!s}
    -totals             Ca Mg Sr Ba C O {totals!s}    
    -molalities         CO2 {molalities!s}
    -isotopes           R(13C)_CO2(aq)  R(18O)_CO2(aq) R(14C)_CO2(aq) R(14C)
                        {isotopes!s}
    -solid_solutions    Calcite
"""

INITIAL_CHEM = """
SOLUTION 1 Soilwater
    temp    {temp:12G}
    units   mmol/kgw
    pH      {ph:12G}    charge
    Ca      {ca:12G}    
    Mg      {mg:12G}
    Sr      {sr:12G}
    Ba      {ba:12G}
    Cl      {cl:12G}    
    C       {c:12G}    # calculate using _get_init_co()
    O(0)    {o:12G}    # calculate using _get_init_co()
    [18O]   {d18O:12G}    # permil
    [13C]   {d13C:12G}    # permil
    [14C]   {R14C:12G}    # pmc
    [44Ca]  {d44Ca:12G}    # permil
    -water  1    # kg
"""

DUMMY_CHEM_CO = """
SOLUTION 0 Determine C and O concentrations
    temp    {temp:12G}
    units   mmol/kgw
    pH      {ph:12G}    charge
    Ca      {ca:12G}    
    Mg      {mg:12G}
    Sr      {sr:12G}
    Ba      {ba:12G} 
    Cl      {cl:12G}    
    C       1    CO2(g)    {co2_si:12G}
    O(0)    1    O2(g)     {o2_si:12G}
    -water  1    # kg
"""

DUMMY_CHEM_D13C = """
SOLUTION 0 Determine solution d13C from gas
    temp    {temp:12G}
    units   mmol/kgw
    pH      {ph:12G}    charge
    Ca      {ca:12G}    
    Mg      {mg:12G}
    Sr      {sr:12G}
    Ba      {ba:12G} 
    Cl      {cl:12G}    
    C       {c:12G}    # calculate using _get_init_co()
    O(0)    {o:12G}    # calculate using _get_init_co()
    [18O]   {d18O:12G}    # permil
    [13C]   {d13C:12G}    # permil
    [14C]   {R14C:12G}    # pmc
    [44Ca]  {d44Ca:12G}   # permil
    -water  1    # kg
"""

WATER_EQUILIBRATE = """
REACTION            # force isotopic equilibration
    H2O(g)   0
"""


BEDROCK_PHASES = """
PHASES
Bedrock
    {!s}
"""

# BEDROCK_DISS 1 uses Calcite saturation (more geochemically realistic)
# for use with quantitative dissolutions
BEDROCK_DISS1 = """
EQUILIBRIUM_PHASES
    Pyrite  0   {moles_pyrite}  dissolve_only
    Calcite 0   Bedrock {moles_bedrock} dissolve_only
        -force_equality


"""

# BEDROCK_DISS 2 uses bedrock dissolution parameters
# for use with incongruent dissolution (i.e. where the name 'Calcite' is used
# for a precipitating reaction, and hence can't be used here).
BEDROCK_DISS2 = """
EQUILIBRIUM_PHASES
    Bedrock 0   {moles_bedrock} dissolve_only
    Pyrite  0   {moles_pyrite}  dissolve_only
        -force_equality


"""

CO2_OUT = """
REACTION
    CO2         {molefrac_co2:12G}
    [13C]O2     {molefrac_13C:12G}
    CO[18O]     {molefrac_18O:12G}
    [14C]O2     {molefrac_14C:12G}
    {gas_moles:12G}
"""

CO2_OUT2 = """
REACTION
    CO2         {molefrac_co2:12G}
    [13C]O2     {molefrac_13C:12G}
    CO[18O]     {molefrac_18O:12G}
    {gas_moles:12G}
"""

CO2_GAS = """
GAS_PHASE {name}
    -fixed_pressure
    -pressure       {pressure:12G}
    -volume         {volume:12G}
    -temperature    {temperature:12G}
    CO2(g)          {pp_co2:12G}
    CO[18O](g)      {pp_C18O:12G}
    C[18O]2(g)      {pp_C18O2:12G}
    [13C]O2(g)      {pp_13C:12G}
    [13C]O[18O](g)  {pp_13C18O:12G}
    [13C][18O]2(g)  {pp_13C18O2:12G}
    [14C]O2(g)      {pp_14C:12G}
    [14C]O[18O](g)  {pp_14C18O:12G}
    [14C][18O]2(g)  {pp_14C18O2:12G}
    O2(g)           {pp_O2:12G}
    O[18O](g)       {pp_O18O:12G}
    Q2(g)           {pp_N2:12G}

"""

CALCITE_SS = """
SOLID_SOLUTION 1 Calcite
    Calcite
    -comp Calcite 0
    -comp CaCO2[18O](s) 0
    -comp CaCO[18O]2(s) 0
    -comp CaC[18O]3(s) 0
    -comp Ca[13C]O3(s) 0
    -comp Ca[13C]O2[18O](s) 0
    -comp Ca[13C]O[18O]2(s) 0
    -comp Ca[13C][18O]3(s) 0
    -comp [44Ca]CO3(s) 0
    -comp [44Ca]CO2[18O](s) 0
    -comp [44Ca]CO[18O]2(s) 0
    -comp [44Ca]C[18O]3(s) 0
    -comp [44Ca][13C]O3(s) 0
    -comp [44Ca][13C]O2[18O](s) 0
    -comp [44Ca][13C]O[18O]2(s) 0
    -comp [44Ca][13C][18O]3(s) 0
    # -comp Ca[14C]O3(s) 0   # radiocarbon may cause non-convergence
    # -comp Ca[14C]O2[18O](s) 0
    # -comp Ca[14C]O[18O]2(s) 0
    # -comp Ca[14C][18O]3(s) 0
    # -comp [44Ca][14C]O3(s) 0
    # -comp [44Ca][14C]O2[18O](s) 0
    # -comp [44Ca][14C]O[18O]2(s) 0
    # -comp [44Ca][14C][18O]3(s) 0            
"""

SS_MG = """   # Mg components:         -
    -comp MgCO3(s) 0
    -comp MgCO2[18O](s) 0
    -comp MgCO[18O]2(s) 0
    -comp MgC[18O]3(s) 0
    -comp Mg[13C]O3(s) 0
    -comp Mg[13C]O2[18O](s) 0
    -comp Mg[13C]O[18O]2(s) 0
    -comp Mg[13C][18O]3(s) 0
    # -comp Mg[14C]O3(s) 0
    # -comp Mg[14C]O2[18O](s) 0
    # -comp Mg[14C]O[18O]2(s) 0
    # -comp Mg[14C][18O]3(s) 0

"""

SS_SR = """   # Sr components:
     -comp SrCO3(s) 0
    -comp SrCO2[18O](s) 0
    -comp SrCO[18O]2(s) 0
    -comp SrC[18O]3(s) 0
    -comp Sr[13C]O3(s) 0
    -comp Sr[13C]O2[18O](s) 0
    -comp Sr[13C]O[18O]2(s) 0
    -comp Sr[13C][18O]3(s) 0
    # -comp Sr[14C]O3(s) 0
    # -comp Sr[14C]O2[18O](s) 0
    # -comp Sr[14C]O[18O]2(s) 0
    # -comp Sr[14C][18O]3(s) 0
"""

SS_BA = """   # Ba components:
     -comp BaCO3(s) 0
    -comp BaCO2[18O](s) 0
    -comp BaCO[18O]2(s) 0
    -comp BaC[18O]3(s) 0
    -comp Ba[13C]O3(s) 0
    -comp Ba[13C]O2[18O](s) 0
    -comp Ba[13C]O[18O]2(s) 0
    -comp Ba[13C][18O]3(s) 0
    # -comp Ba[14C]O3(s) 0
    # -comp Ba[14C]O2[18O](s) 0
    # -comp Ba[14C]O[18O]2(s) 0
    # -comp Ba[14C][18O]3(s) 0    
"""