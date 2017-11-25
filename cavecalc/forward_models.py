"""Allows multiple models to be run in sequence and handles model IO. For most
use-cases, cavecalc Models should be run using the ForwardModels class.

Classes defined here:
    ForwardModels
    
Functions defined here:
    run_a_model
    run_async
    run_linear
"""

# built-in modules
import os
import multiprocessing as mp
import pickle
# my modules
import cavecalc.caves as ccv
import cavecalc.util as ccu
from cavecalc.setter import SettingsObject, SettingsMaker
from copy import deepcopy

def run_a_model(SO):
    """Run a single cavecalc model.
    
    Args:
        SO (SettingsObject): A SettingsObject describing model input.
    Returns:
        (r, id) where r is the model results dict and id is the 'id' parameter
        in settings. If the model raises an error, (None, id) is returned.
    """
    
    i = SO.id
    sim = ccv.Simulator(SO.dict(), i)
    r = sim.run()
    print("Model %i complete." % i)
    return (r, i)
        
def run_async(SO_list):
    """Run multiple models in parallel.
    
    WARNING: Parallel processing may cause errors in some calculations. Use
    only when absolutely necessary. run_linear() is recommended instead.
    
    Args:
        SO_list: A list SettingsObjects, e.g. from a SettingsMaker
    Returns:
        A list of (r, id) tuples. r is the model results dict and id is the 
        'id' parameter in the settings dict.
    """
    
    s = SO_list
    results = {}
    pool = mp.Pool()
    results = [pool.apply_async(run_a_model, (e.dict(), e.id)) for e in s]
    pool.close()
    pool.join()
    r = [res.get() for res in results]
    r2 = sorted(r, key=lambda tup: tup[1]) # sort by model id
    return r2
    
def run_linear(SO_list):
    """Runs multiple models in sequence.
    
    Args:
        SO_list: A list SettingsObjects, e.g. from a SettingsMaker
    Returns:
        A list of (r, id) tuples. r is the model results dict and id is the 
        'id' parameter in the settings dict.
    """    
    
    return [run_a_model(e) for e in SO_list]

class ForwardModels(object):
    """Runs Cavecalc models.
    
    Handles generation and checking of settings suites, and saving of bundled
    model output. This class provides a flexible interface to cavecalc.caves 
    (the module which actually runs the models), and is the preferred method 
    of running models.
    """
    
    def __init__(self, settings=None, output_dir=None):
        """Initialise the object and process settings.
        
        Args:
            settings (SO): Model input parameters. See defaults.txt for 
                             options.
            settings_file: If no settings arg is given, settings may be loaded
                           from a file.
            output_dir: (Optional) Specify path to directory for saving files.
                        By default files are saved to the current directory.
        """
        
        self.done_input = []
        self.done_results = []
        
        # Get settings objects
        if settings:    self.input = SettingsMaker(**settings).settings()
        else:           self.input = SettingsMaker().settings()
            
        if output_dir:  
            if not os.path.isdir(output_dir): 
                os.mkdir(output_dir)
            self.output_dir = output_dir
        else: self.output_dir=os.getcwd()

    def _check_previous_saves(self, interactive=False, use_by_default=True):
        """Checks output directory for existing output and prompts the user
        to decide whether they want to use it or not."""
        
        
        def dict_find(SO, list_of_SOs):
        
            for i,s in enumerate(list_of_SOs):
            
                # check keys are equal
                if SO.dict().keys() != s.dict().keys():
                    continue
                    
                # check all values are equal
                equal = True
                for k,v in SO.dict().items():
                    if v != s.dict()[k]:
                        equal = False
                if equal:
                    return i
            return None 
        
        try:
            with open(os.path.join(self.output_dir, 'settings.pkl'), 'rb') as f:
                prev_input = pickle.load(f)
                
            with open(os.path.join(self.output_dir, 'results.pkl'), 'rb') as f:
                prev_results = pickle.load(f)
        except FileNotFoundError:
            return

        if (prev_input is None) or (prev_results is None):
            return
        
        new_input = []
        done_input = []
        done_results = []
        
        # for i,prev in enumerate(prev_inputs):
        for s in self.input:
            i = dict_find(s, prev_input)
            if i is None:
                new_input.append(s)
            else:
                done_input.append(prev_input[i])
                done_results.append(prev_results[i])

        print("Previous model output detected for selected input settings.")
        
        if interactive == True:
            a = ''
            while a == '':
                print( "%i out of %i models appear to be repeated." % 
                       (len(done_input), len(self.input)))
                a = input("Re-use old output for these models (y/n)? ")
                if a.capitalize() == 'Y':
                    reuse = True
                elif a.capitalize() == 'N':
                    reuse = False
                else:
                    a = ''
                
        else:
            print("%i out of %i models are repeats." % 
                       (len(done_input), len(self.input)))
            if use_by_default:
                print("Re-using old calculations where available.")
            else:
                print("Re-running all models")
            reuse = use_by_default
        
        if reuse:
            self.input = new_input
            self.done_input = done_input
            self.done_results = done_results
    
    def run_models(self, mode='Serial', **kwargs):
        """Run models for all parameter sets loaded into the object. Output is
        addded to self.results as a list of (r, id) tuples. See run_linear and
        run_async for method. Parallel mode is not currently functional (on
        Windows at least).
        
        Args:
            mode: 'Serial' (default) or 'Parallel'. Running models in parallel
                  is faster but not recommended - results may be inaccurate, as
                  Iphreeqc is not necessarily thread-safe in this 
                  implementation.
        Returns:
            Nothing. Model results list is assigned to self.results. List 
            indices in self.results correspond to indices in self.self.input.
        """

        self._check_previous_saves(**kwargs)
        ret_dir = os.getcwd()
        print("Models to run:\t%s" % len(self.input))
        try:
            os.chdir(self.output_dir)
            if mode.capitalize() == 'Serial': # Run models in order
                results = run_linear(self.input)
            elif mode.capitalize() == 'Parallel': # Run models multithreaded
                print("Warning: Parallel operation accuracy not guaranteed.")
                results = run_async(self.input)
            else:
                raise ValueError("Mode not recognised. Use serial / parallel")
            self.results = results
            
            
            # add re-used output, if any
            self.results.extend(self.done_results)
            self.input.extend(self.done_input)
            
        finally:
            os.chdir(ret_dir)
            
    def save(self):
        """Save results and settings data to .pkl files. 
        
        The resulting files (results.pkl and settings.pkl) may be read using 
        the cavecalc.analyse module.
        
        results.pkl contains a list of dicts. Each dict contains the output of
        a single model. settings.pkl contains a similar list of 
        SettingsObjects, one for each model run.
        """
        
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        with open(os.path.join(self.output_dir, 'settings.pkl'), 'wb') as f:
            pickle.dump(self.input, f)
            
        with open(os.path.join(self.output_dir, 'results.pkl'), 'wb') as f:
            pickle.dump(self.results, f)
        
    def _debug(self, i):
        """Runs the specified model and saves the pq_input_log file for 
        debugging. The input log file may be run directly by PHREEQC.
        
        Args:
            i: Index position in self.input to re-run in debug mode.
        """
        s_debug = self.input[i]
        s_debug.set(phreeqc_log_file=True)
        
        dbug = ccv.Simulator(s_debug, id=i)
        dbug.run()