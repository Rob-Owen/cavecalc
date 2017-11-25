"""Contains the Evaluate object, which has methods to load data from Cavecalc
.pkl output files, display data and write it to other file formats.

Classes defined here:
    Evaluate
"""

import pickle
import os
import copy
import matplotlib
from sys import platform
if platform != 'win32':
    matplotlib.use('TkAgg') # necessary for mac
from matplotlib import pyplot as plt
import numpy as np
import cavecalc.util as ccu
import scipy.io as sio
import seaborn as sns

class Evaluate(object):
    """Processes of Cavecalc model output.
    
    Evaluate contains methods to load, process, format and display Cavecalc
    model data saved in .pkl files. It allows model output to be saved to .csv
    and .mat files for further processing.
    
    Evaluate objects may be used directly (see examples), and are also used
    under the hood by the output GUI.
    
    Methods:
        get_settings_report - Get a dict summarising model settings
        load_data           - Load .pkl files from a directory
        save_csvs           - Save all loaded model output to .csv files
        save_all_mat        - Save all loaded model output to a .mat file
        filter_out_noprecip - Filter out all model steps that don't
                              precipitate calcite.
        filter_by_index     - Filter model steps in/out based on a model step 
                              index.
        filter_by_results   - Filter model steps in/out based on any output 
                              parameter.
        filter_by_settings  - Filter whole models in/out based in any input 
                              parameter.
        plot_models         - Plot all loaded models.
        plot_points         - Plot selected steps from all loaded models.
    """
    
    def __init__(self):
        """Initialise an Evaluate object.
        
        After initialisation, load_data must be called to read data from .pkl
        files.
        """
        
        self._models = []
        self._settings = []

    @property
    def model_settings(self):
        """Return a list of dicts containing model settings.
        
        Generates a list of dicts from from all loaded SettingsObjects by
        calling their .dict() method.
        
        Returns:
            A list of settings dicts.
        """
        
        if not self._settings:
            raise ValueError("Object %r has no models loaded." % self)
        
        o = [s.dict() for s in self._settings]
        for d in o:
            try:
                d.pop('id')
            except KeyError:
                pass
            
        if o:
            return copy.deepcopy(o)
        else:
            raise ValueError("Object %r has no models loaded." % self)
    
    @property
    def model_results(self):
        """Return a list of dicts containing model output.
        
        Returns:
            A list of results dicts.
        """
        
        if self._models:
            return copy.deepcopy(self._models)
        else:
            raise ValueError("Object %r has no models loaded." % self)
    
    def get_settings_report(self):
        """Get a summary of the range of model settings.

        Returns:
        A dict of model settings, with one entry for each unique value detected.
        """

        d = self.model_settings[0]
        o = dict.fromkeys(d.keys(),[])
        for s in self.model_settings:
            for k, v in s.items():
                if v not in o[k]:
                    o[k].append(v)
        try:
            o.pop('id')
        except KeyError:
            pass
        return o

    def load_data(self, *args):
        """Load .pkl data into the Evaluate object for processing.
        
        Data is loaded from a directory. The directory must contain
        settings.pkl and results.pkl. load_data may be called multiple times
        to merge different model suites for comparison.
        
        Args:
            *args: The directories to load data from.
        """        

        ret_dir = os.getcwd()
        
        if len(args) == 0:
            args = (ret_dir,)
            
        for d in args:
            os.chdir(d)
            try:
                print("Attempting to load data from %s..." % d, end="")
                with open('settings.pkl', 'rb') as f:
                    self._settings.extend(pickle.load(f))

                with open('results.pkl', 'rb') as f:
                    r = pickle.load(f)
                    self._models.extend([a for (a,b) in r])
                print(" Done")
            finally:
                os.chdir(ret_dir)

    def save_csvs(self, directory=None):
        """Save model output to .csv files.
        
        One file is saved for each model loaded. Note that only model output
        can be saved to csv, not the settings used.
        
        csv files may be easily read in any spreadsheet program.
        
        Args:
            directory (str): The directory to save output to.
        """
        
        if not directory:
            directory = os.getcwd()
        
        for (id, model) in enumerate(self._models):
            f = os.path.join(directory, "out_%i.csv" % id)
            ccu.save_csv(model, os.path.join(f))

    def save_all_mat(self, file):
        """Save all loaded data to a .mat file.
        
        Data is saved as two matlab structs, reflecting the data structures
        inside settings.pkl and results.pkl respectively.
        
        Args:
            file: Filename to save to (.mat will be auto-appended)
        """

        s = dict()
        
        for i, SO in enumerate(self._settings):            # for each model
            set = SO.dict()         # settings dict
            res = self._models[i]    # results dict
            
            # remove any 'None' values from output (savemat can't handle them)
            # replace with -999 value, like PHREEQC
            n_res = dict()
            for k,v in res.items():
                nv = []
                for e in v:
                    if e is None:
                        nv.append(-999)
                    else:
                        nv.append(e)
                n_res[k] = nv
        
            o = {k:(v if type(v) is list else [v]) for k,v in set.items()}
            
            a = ccu.numpify(o)                          # settings
            b = ccu.numpify(n_res)                      # output
            
            c = { 'settings' : a,
                  'results'  : b   }
            
            name = "m" + str(i)
            s[name] = c

        sio.savemat(file, s)

    def filter_out_noprecip(self):
        """Returns a filtered copy of the Evalulate object.
        
        Models are filtered out of they do not include any precipitation
        reactions. This is useful for 'eq' mode analyses to remove 
        non-precipitating solutions.
        """
        
        A = Evaluate()
        
        for i,m in enumerate(self._models):
            a = False
            for s in m['step_desc']:
                if 'precip' in s:
                    a = True
                    break
            
            if a:
                A._models.append(copy.deepcopy(m))
                A._settings.append(copy.deepcopy(self._settings[i]))
        
        return A
    
    def filter_by_index(self, ind, n=False):
        """Return a filtered copy of the Evaluate object.
        
        Filter the model output data contained in the object. Data is filtered
        based on list index position - this corresponds to the calculation step
        in the model. This method is useful for subsetting data in preparation 
        for plotting. It works similarly to filter_by_results.
        
        Example:
            e = Evaluate()
            e.load_data('./my_data/')
            f = e.filter_by_index(-1) # extracts the final dripwater chemistry
           
        Args:
            ind: An integer index to filter by. This corresponds to a model
                step number. E.g. index 0 is the first PHREEQC calculation
                (initial water chemistry), index 1 is the bedrock dissolution
                product, index -1 is the final solution chemistry.
            n: Optional boolean argument. If True, the filter is inverted.
                Default False.
                
        Returns:
            A modified copy of the object. The copy only contains model output
            that meet the filter criteria.
        """
        A = Evaluate()
        A._settings = copy.deepcopy(self._settings)
        
        for m in self._models:
            if ind < 0:
                explicitIndex = len(m['step_desc']) + ind
            else:
                explicitIndex = ind

            if n is False:
                fil = {k : [a[explicitIndex]] for k,a in m.items()}
            else:
                fil = {k : [v for i,v in enumerate(a) if i != explicitIndex] for k,a in m.items()}
            A._models.append(fil)
            
        
        rem = []
        for i,r in enumerate(A._models):
            if max([len(v) for k,v in r.items()]) == 0:
                rem.append(i)
                
        [A._models.pop(i) for i in rem]
        [A._settings.pop(i) for i in rem]
        
        return copy.deepcopy(A)
            
    def filter_by_results(self, key, value, n=False):
        """Return a filtered copy of the Evaluate object. 
        
        Filter the model output data contained in the object. Data is filtered
        based on a key, value combination. This method is useful for 
        subsetting data in preparation for plotting. It works similarly to 
        filter_by_index.
        
        Example:
            e = Evaluate()
            e.load_data('./my_data/')
            f = e.filter_by_settings('step_desc', 'degas')
            # f includes only data from the degassing steps
        
        Args:
            key: Key in model output dicts to filter by.
            value: Value to filter 'key' by. Accepts substrings for step_desc.
            n: Optional boolean argument. If True, the filter is inverted.
                Default False.
                
        Returns:
            A filtered copy of the Evaluate object.
        """

        A = Evaluate()
        A._models = []
        A._settings = self._settings

        # filter object
        for i, m in enumerate(self._models):
            fil = {}
            a = m[key]
            for j, v in m.items():
                if len(v) == len(a):
                    if n:
                        fil[j] = [v[k] for k in range(len(v)) if value not in a[k]]
                    else:
                        fil[j] = [v[k] for k in range(len(v)) if value in a[k]]
                else:
                    fil[j] = v
            A._models.append(fil)
        return copy.deepcopy(A)

    def filter_by_settings(self, setting, value, n=False):
        """Return a filtered copy of the Evaluate object.
        
        The returned Evaluate object contains a subset of the models in the
        original. Models are filtered based on the settings, value combination
        provided. Models that meet the critera have their data included in the
        copy.
        
        Args:
            setting (str): model parameter to filter by (e.g. 'gas_volume')
            value: value of 'setting' to include (e.g. 20).
            n: Optional boolean argument. If True, the filter is inverted.
                Default False.
        Returns:
            A filtered copy of the Evaluate object.
        """

        A = Evaluate()
        A._models = []
        A._settings = []

        for i, b in enumerate(self._settings):
            d = b.dict()
            if n:
                if isinstance(value, str):
                    if value not in d[setting]:
                        A._models.append(self._models[i])
                        A._settings.append(b)
                else:
                    if d[setting] != value:
                        A._models.append(self._models[i])
                        A._settings.append(b)
            else:
                if isinstance(value, str):
                    if value in d[setting]:
                        A._models.append(self._models[i])
                        A._settings.append(b)
                else:
                    if d[setting] == value:
                        A._models.append(self._models[i])
                        A._settings.append(b)

        return copy.deepcopy(A)

    def plot_models(self, *args, x_key=None, y_key=None, 
                    label_with=None, ax=None, **kwargs):
        """Plot Model results with one series per model.
        
        Creates a simple matplotlib figure. Useful, for example, to quickly
        display the degassing evolution of a suite of models. May be combined
        with filter_by_settings, filter_by_results or filter_by_index to
        include / exclude certain parts of the dataset.
        
        Args:
            *args: Optional formatting parameters passed to pyplot.plot()
            x_key: Model output to plot on x-axis
            y_key: Model output to plot on y-axis
            label_with (optional): Model input parameter to annotate series
            ax (optional): Add data to a pre-existing matplotlib axis
            **kwargs (optional): kwargs to be passed to pyplot.plot()
        Returns:
            Axes object.
        """
        
        sns.set_style('darkgrid')
        if not ax:
            fig, ax = plt.subplots()
            ax.set_ylabel(y_key)
            ax.set_xlabel(x_key)
        for i, m in enumerate(self._models):
            if label_with:
                s = self._settings[i].get(label_with)
                a = "%s: %s" % (label_with, s)
                ax.plot(m[x_key], m[y_key], label = a, *args, **kwargs)
                ax.legend(prop={'size':6})
            else:
                ax.plot(m[x_key], m[y_key], *args, **kwargs)
            

        return ax

    def plot_points(self, *args, x_key=None, y_key=None, plot_index=1, 
                    label_with=None, ax=None, **kwargs):
        """Plot Model results for a point-by-point inter-model comparison.
        
        Useful, for example, to show different bedrock dissolution products
        across a suite of models.
        
        Args:
            x_key: Model output or setting parameter to plot on x-axis
            y_key: Model output to plot on y-axis
            *args (optional): Formatting parameters passed to pyplot.plot()
            plot_index: Which point to plot. e.g. 0 (initial water), 1 (bedrock
            dissolution product), -1 (fully degassed solution)
            label_with (optional): Model input parameter to label points with
            ax (optional): Add data to a pre-existing plot
            **kwargs (optional): kwargs to be passed to pyplot.plot()
        Returns:
            Axes object.
        """
        sns.set_style('darkgrid')
        
        x_vals = []
        y_vals = []
        labels = []

        # look for x_key in results
        if x_key in list(self._models[0].keys()):
            for i, m in enumerate(self._models):
                try:
                    x_vals.append(m[x_key][plot_index])
                except IndexError:
                    pass

        # otherwise, find it in settings
        else:
            for i, s in enumerate(self._settings):
                x_vals.append(s.dict()[x_key])

        for i, m in enumerate(self._models):
            if label_with:
                s = self._settings[i].dict()
            try:
                y_vals.append(m[y_key][plot_index])
                if label_with:
                    labels.append(s[label_with])
            except IndexError:
                pass
        
        if not ax:
            fig, ax = plt.subplots()
            ax.set_ylabel(y_key)
            ax.set_xlabel(x_key)
        ax.plot(x_vals, y_vals, *args, **kwargs)
        if label_with:
            for lab, x, y in zip(labels, x_vals, y_vals):
                ax.annotate('%s=%s' % (label_with, lab),
                            xy=(x, y), fontsize=8)
        return ax
