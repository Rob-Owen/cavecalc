"""Run a suite of models investigating the effect of temperature on the bedrock
dissolution step under open system conditions. Do it for both Calcite and
Dolomite
"""

from cavecalc.forward_models import ForwardModels
import cavecalc.analyse as cca

#
## STEP 1: Define the non-default settings and run the model.
# use diss_only mode because we don't care about the degassing phase.

s =  {  'bedrock_mineral' :      ['Calcite', 'Dolomite'],
        'temperature'     :      [20, 22, 24, 26, 28, 30],
        'gas_volume'    :        20,
        'init_pCO2' :            15000,
        'init_d13C' :           -15,
        'kinetics_mode' :       'diss_only',
        'cave_pCO2' :           500
     }

dir = './ex3/' # directory to save model output
p = ForwardModels(  settings =      s,
                    output_dir  =   dir )   # initialise the model
p.run_models()                              # run the model
p.save()                                    # save the model output

#
## STEP 2: Load model output and do stuff with it
#
e = cca.Evaluate()      # initialise an 'Evaluate' object
e.load_data( dir )      # load data from dir

#
## STEP 3: Plot pH vs Ca concentration.
#

# plot_points plots a single point from each model
e1 = e.filter_by_settings('bedrock_mineral', 'Calcite')
e2 = e.filter_by_settings('bedrock_mineral', 'Dolomite')

e1.plot_points( x_key='Ca(mol/kgw)', y_key='pH', 
               label_with = 'temperature',
               plot_index = 1 )
               
a1 = e1.plot_points( x_key='temperature', y_key='Ca(mol/kgw)', 
               label_with = 'bedrock_mineral',
               plot_index = 1 )
               
e2.plot_points( x_key='temperature', y_key='Ca(mol/kgw)', 
               label_with = 'bedrock_mineral',
               plot_index = 1, ax=a1 )
cca.plt.show()
