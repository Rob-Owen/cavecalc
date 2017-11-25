"""
Runs Cavecalc with default settings and:
    1) Save a csv file with the results
    2) Make a plot
"""

from cavecalc.forward_models import ForwardModels
import cavecalc.analyse as cca

#
## STEP 1: Run the model
#
dir = './ex1/' # directory to save model output
p = ForwardModels( output_dir = dir )    # initialise the model with default inputs
p.run_models()                              # run the model
p.save()                                    # save the model output

#
## STEP 2: Load model output and do stuff with it
#
e = cca.Evaluate()      # initialise an 'Evaluate' object
e.load_data( dir )      # load data from dir
e.save_csvs( dir )      # save the model output to a .csv file

#
## STEP 3: Make a plot of pH against Ca concentration. Use red crosses.
#
e2 = e.filter_by_index(0, n=True) # filter out the first model step
e2.plot_models('-xr', x_key='Ca(mol/kgw)', y_key='d13C')
cca.plt.show()
