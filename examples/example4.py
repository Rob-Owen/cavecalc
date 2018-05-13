"""A more advanced API example - this is not possible using the GUI.

Here we undergo a degassing loop similar to running Cavecalc in 'eq' mode, 
but also allow some amount of equilibration with the bedrock on each model 
step. 3 models are run, each with a different amount of equilibration allowed.

In previous examples we use the ForwardModels object to handle the running of
the model. Here we skip that step and directly access the objects in caves.py 
to allow us more flexibility.

Consult the documentation for caves.py for detailed info on all methods and
objects used.
"""

from cavecalc.caves import *
from cavecalc.setter import SettingsObject
import cavecalc.util as ccu
import matplotlib.pyplot as plt

s =  {  'bedrock_mineral' :      'Calcite',
        'temperature'     :      32,
        'gas_volume'    :        20,
        'init_pCO2' :            15000,
        'init_d13C' :           -15,
        'phreeqc_log_file' :    True,
        'cave_pCO2' :           500
     }

# we run three models, each with a different amount of bedrock equilibration 
for m_bedrock in [0, 1e-5, 1e-3]: # three values for m (moles of bedrock)
    
    # A SettingsObject fills in s with default values for all other parameters
    # and checks our specified inputs are of valid types.
    SO = SettingsObject(**s)
    
    # sim - a Simulator object - is the core of the model. It handles all IO
    # from IPhreeqc. We initialise one by providing it with a SettingsObject.
    sim = Simulator(SO)

    # initialise bedrock and water objects. These are associated with the
    # simulator 'sim' and get their initial settings from it.
    bedrock = Carbonate(sim)
    water = Solution(sim)
    
    # Here we do the first reaction. water.set_soil_chem() returns PHREEQC
    # input to define the initial soil chemistry. This is added to the Iphreeqc
    # buffer, to be executed later by sim.ipq_exec().
    sim.ipq_buffer( *water.set_soil_chem() )
    sim.ipq_exec()    #    Calculate initial soilwater chemistry

    # Next reaction: dissolve bedrock in the soilwater.
    sim.ipq_buffer(*water.open_diss_rxn( bedrock ))
    sim.ipq_exec()    #    Calculate bedrock dissolution product

    # Now do an interative degassing loop until 'cave_pCO2' is reached.
    while water.will_it_degas( sim.settings['cave_pCO2'] ):
        # on each loop, add 3 reactions to the buffer.
        # water.wri_rxn is the one not used in the standard 'eq' mode
        sim.ipq_buffer([ water.wri_rxn(bedrock, m_bedrock=m_bedrock),
                         water.kinetic_degas_rxn(),
                         water.precipitate_rxn() 
                       ])
        sim.ipq_exec()

    ccu.PostProcessor(sim) # perform offline calculations and processing

    #### look at output
    x = 'f_ca'
    y = 'd13C'
    plt.figure()
    plt.plot(sim.output[x][1:], sim.output[y][1:])
    plt.xlabel(x)
    plt.ylabel(y)

plt.show()