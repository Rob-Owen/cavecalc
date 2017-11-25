"""
Run two models, one with dolomite bedrock and one with calcite bedrock. All
other parameters are set to default values.

Also save phreeqc log file recording how the model ran.
"""

from cavecalc.forward_models import ForwardModels
import cavecalc.analyse as cca

#
## STEP 1: Define the non-default settings and run the model
#

s =  {  'bedrock_mineral' :      ['Calcite', 'Dolomite'],
        'phreeqc_log_file':      True,
        'cave_pCO2' :            1000
     }

dir = './ex2/' # directory to save model output
p = ForwardModels(  settings =      s,
                    output_dir  =   dir )   # initialise the model
p.run_models()                              # run the model
p.save()                                    # save the model output

#
## STEP 2: Load model output and save new data formats
#
e = cca.Evaluate()      # initialise an 'Evaluate' object
e.load_data( dir )      # load data from dir
e.save_all_mat( dir )   # save data to a .mat file
e.save_csvs( dir )      # save each model to a .csv file

#
## STEP 3: Plot pH vs Ca concentration.
#
e2 = e.filter_by_index(0, n=True) # filter out the first model step
e2.plot_models(x_key='Ca(mol/kgw)', y_key='pH', label_with = 'bedrock_mineral')
cca.plt.show()
