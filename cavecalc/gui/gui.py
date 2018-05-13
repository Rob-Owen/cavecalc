"""Codes for the Cavecalc GUIs.

The GUIs provide a simple user-friendly interface for Cavecalc. They do not
expose the full model functionality but allow most common calculations to be
performed.

GUIs may be accessed by running the scripts:
    cc_input_gui.py
    cc_output_gui.py
"""

import os
import copy
from collections import OrderedDict
import operator
from numpy import linspace
import matplotlib
from sys import platform
if platform != 'win32':
    matplotlib.use('TkAgg') # necessary for mac
from matplotlib import pyplot as plt
from tkinter import *
from tkinter import filedialog
from cavecalc.analyse import Evaluate
import cavecalc.data.types_and_limits
import cavecalc.gui.mapping  
import cavecalc.gui.layout
from cavecalc.setter import SettingsMaker, NameSwitcher, SettingsObject

# settings options hidden from plotting menus (useless clutter)
HIDE_OPS = True
HIDDEN_OPTS = [
    'totals',  'molalities',   'isotopes',     'out_dir', 
    'phreeqc_log_file',        'phreeqc_log_file_name',
    'database'] # options not available for plotting

ns = NameSwitcher()

def od(dict):
    """Converts a dict to an ordered dict, sorted by key."""
    
    return OrderedDict(sorted(dict.items()))

def py2tk(dict):
    """Converts a dict to Tkinter types.
    
    Convert dictionary entries from doubles and strings to StringVar for use 
    with tkinter. Returns a modified copy.
    
    Args:
        dict: A dict with entries that are simple data types.
    Returns:
        A modified dict.
    """
    
    out = copy.copy(dict)
    types = vars(cavecalc.data.types_and_limits).copy()
    
    for k in dict.keys():
        if type(types[k]) is bool:
            out[k] = BooleanVar()
            out[k].set(False)
        elif dict[k] is not None:
            out[k] = StringVar()
            out[k].set(dict[k])
        else:
            out[k] = None
    return out
    
def _parse_value_input(string, allow=[]):
    """
    Parses leftmost panel input to detect ranges of values or single values.
    Returns either a double or a list of doubles.
    """
    a = string
    
    # remove brackets
    rm = ['(', ')', '[', ']']
    for s in rm:
        if s not in allow:
            a = a.replace(s,'')
    
    # replace commas with space
    a = a.replace(',',' ')
    
    #split string on spaces
    a = a.split(' ')
    
    # attempt conversion to float (most data types are numeric)
    try:
        b = [float(v) for v in a if v != '']
    except ValueError:
        b = [str(v) for v in a if v!= '']
    
    # remove list structure from single entries
    if len(b) == 1:
        b = b[0]
    
    return b
    
def tk2py(dict, parse=False):
    """Inverse of py2tk.
    
    Convert a dict of tkinter StringVar types to a a dict understandable by
    the cavecalc setter module.
    
    Args:
        dict: A dict full of GUI inputs (e.g. StringVar types)
        parse (optional): If True, process numeric input to get a list of
            values if possible. Default False.
    Returns:
        A dictionary of booleans, strings, lists and floats.
    """    
    
    a = dict.copy()
    for k in dict.keys():
        if dict[k] is None:
            a[k] = None
        elif type(dict[k].get()) is str:
            if parse:
                a[k] = _parse_value_input(dict[k].get())
            else:
                a[k] = dict[k].get()
        else:
            a[k] = dict[k].get()
        
    b = {k:v for k,v in a.items() if v is not None}
    return b

def gplot(  x_values, y_values, x_label, y_label, 
           label_vals, label_name ):
    """Plot data in a new window.
    
    This plotting function is a simplified version of the plotting code used in
    the analyse module for use with the GUI. If x_values contains lists of
    length 1, a single series is plotted, connecting data points from all
    models together. If length > 1 then each model is plotted as a separate 
    series.
    
    If more advanced plotting is required... use something else.
    
    Args:
        x_values: A list of lists containing the model output to be plotted on
            the x-axis.
        y_values: A list of lists containing the model output to be plotted on
            the y-axis.
        x_label (str): x-axis label.
        y_label (Str): y-axis label.
        label_vals = A list of values to label the series with.
        label_name (str) = The name of the data in label_vals.
    
    """

    fig, ax = plt.subplots()
    
    # if plotting... a single point from each model
    if all(len(x) == 1 for x in x_values):
        xs = [x[0] for x in x_values]
        ys = [y[0] for y in y_values]
        
        ax.plot(xs, ys, 'x--')
        
        if label_name:
            for label, x, y in zip(label_vals, xs, ys):
                ax.annotate(    "%s=%s" % (label_name, label), 
                                xy = (x,y), fontsize=8)
        plt.ylabel(y_label)
        plt.xlabel(x_label)
        
        plt.show(block=False)
     
    # else plot each model as it's own series
    else:
        for (i, xs) in enumerate(x_values):
            ys = y_values[i]
            if label_name:
                l_str = "%s: %s" % (label_name, label_vals[i])
                ax.plot(xs, ys, label = l_str)
                ax.legend(prop={'size':8})
            else:
                ax.plot(xs, ys)
            
        plt.ylabel(y_label)
        plt.xlabel(x_label)
        
        plt.show(block=False)
    
class FileFindWidget(Frame):
    """A Tkinter Widget for opening a file browser window."""
    
    def __init__(self, master=None, value=None, mode=None):
        super().__init__(master)
        self.master = master
        self.value = value
        
        self.entry = Entry(self, textvariable=value).grid(row=0, column=0)
        
        if mode.capitalize() == 'Load':
            self.button = Button(self, text='browse', command= self._openfilename)
        elif mode.capitalize() == 'Save':
            self.button = Button(self, text='browse', command= self._saveasfilename)
        elif mode.capitalize() == 'Dir':
            self.button = Button(self, text='browse', command= self._getdirectory)
        else:
            raise ValueError("Mode %s not recognised. Use save or load." % mode)
        self.button.grid(row=0, column=1)
        
    def _openfilename(self, event=None):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.value.set(file_path)
        
    def _saveasfilename(self, event=None):
        file_path = filedialog.asksaveasfilename()
        if file_path:
            self.value.set(file_path)
            
    def _getdirectory(self, event=None):
        dir_path = filedialog.askdirectory()
        if dir_path:
            self.value.set(dir_path)
             
class InputsRangeWidget(Frame):
    """A compound widget for inputting a range of numeric values.

    The InputsRangeWidget contains an Entry field for inputting a single (or
    multiple values) manually. It also has a 'Set Range' button that opens
    a new window for generating an equally-spaced sequence.
    """
    
    def __init__(self, master=None, value=None):
        super().__init__(master)
        self.master = master
        self.value = value
        self.min = DoubleVar()
        self.max = DoubleVar()
        self.steps = DoubleVar()
        
        self.entry = Entry(self, textvariable=self.value)
        self.entry.grid(row=0, column=0)
        self.button = Button(self, text='*', command=self.select_range)
        self.button.grid(row=0, column=1)
    
    def select_range(self):
        """Creates new window to input range information."""
        
        def get_range():
            vals = linspace( self.min.get(), self.max.get(), 
                             num = self.steps.get())
            self.value.set(vals.tolist())
            top.destroy()
            
        top = Toplevel(self)
        top.title("Select Range") # NOT WORKING / NOT SHOWING
        Label(top, text="Select Range").grid(row=0, columnspan=2)
        Label(top, text="Min").grid(row=1)
        Label(top, text="Max").grid(row=2)
        Label(top, text="Steps").grid(row=3)
        e1 = Entry(top, textvariable = self.min).grid(row=1,column=1)
        e2 = Entry(top, textvariable = self.max).grid(row=2,column=1)
        e3 = Entry(top, textvariable = self.steps).grid(row=3,column=1)
        Button(top, text="OK", command=get_range).grid(row=4)
        Button(top, text="Cancel", command=top.destroy).grid(row=4, column=1)

class OptsWidget(Frame):
    def __init__(self, master, parameter, value, options, row):
        super().__init__(master)
        self.r = row
        self.master = master
        self.parameter = parameter
        self.value = value
        self.options = options
        
        assert type(options) is list
            
        self.l = Label(master, text=self.parameter)
        self.l.grid(row=self.r, column=0, sticky=W)
        self.o = OptionMenu(master, self.value, *self.options)
        self.o.grid(row=self.r, column=1)
        
class InputsWidget(Frame):
    def __init__(self, master, parameter, value, row):
        super().__init__(master)
        self.r = row
        self.master = master
        self.parameter = parameter
        self.value = value
        
        self.l = Label(master, text=parameter)
        self.l.grid(row=self.r, column=0, sticky=W)
        self.o = Entry(master, textvariable=value)
        self.o.grid(row=self.r, column=1)
             
class CCInputGUI(object):
    """The Cavecalc input GUI window."""
    
    def __init__(self, master):
        """Initialise - Open the window."""
        
        self.master = master
        self.master.title('Cavecalc Model Input GUI')
        self._load_defaults()
        self.construct_inputs()
        
    def _loop_gen(self, layout_numbers):
        ln = layout_numbers
        if type(ln) is int:
            ln = [ln]
        elif type(ln) is not list:
            raise TypeError("layout_numbers must be int or list of ints.")
                        
        g = [(k, self.layout[k][1]) for k in self.settings.keys() 
             if self.layout[k][0] in ln]
        
        return sorted(g, key=lambda tup: ns(tup[0]))
        
    def get_ln(self, key):
        """Gets the layout index number of 'key'"""
        
        try:
            return self.layout[key][1]
        except KeyError:
            return self.layout[ns(key)][1]
        
    def _load_defaults(self):
    
        self.d = SettingsObject()

        self.units = vars(cavecalc.data.types_and_limits).copy()
        self.layout = vars(cavecalc.gui.layout).copy()
        
        settings = self.d.dict()
        self.settings = py2tk(settings)
        
    def construct_inputs(self):
        """Frame 1 contains the left-hand panel of the input GUI."""
       
        def add_things_to_frame(frame, layout_number, header_text, i):
            l = Label(frame, text=header_text, font="-size 13")
            l.grid(row=i,columnspan=2, sticky=SW, pady=3)
            i += 1
            for a, b in self._loop_gen(layout_number):
                if b == 'A': # text with range
                    Label(frame, text=ns(a)).grid(row=i, sticky=W)
                    x = InputsRangeWidget(frame, self.settings[a])
                    x.grid(row=i, column=1, sticky=W)
                elif b == 'B': # text without range
                    Label(frame, text=ns(a)).grid(row=i, sticky=W)
                    x = Entry(frame, textvariable=self.settings[a], width=25)
                    x.grid(row=i, column=1, columnspan=2, sticky=W)    
                elif b == 'C': # options menu
                    x = OptsWidget( frame, ns(a), self.settings[a], 
                                    self.units[a], row=i )
                    x.grid(row=i, column=0, sticky=W)
                elif b == 'D': # check button
                    l = Label(frame, text=ns(a)).grid(row=i, sticky=W)
                    r = Checkbutton( frame, variable=self.settings[a],
                                     onvalue=True, offvalue=False )
                    r.grid(row=i, column=1)
                elif b == 'E': # load button
                    l = Label(frame, text=ns(a)).grid(row=i, sticky=W)
                    f = FileFindWidget( frame, value=self.settings[a], 
                                        mode='load')
                    f.grid(row=i, column=1)
                elif b == 'F': # save button
                    l = Label(frame, text=ns(a)).grid(row=i, sticky=W)
                    f = FileFindWidget( frame, value=self.settings[a], 
                                        mode='dir')
                    f.grid(row=i, column=1)
                i += 1
            return i
        
        px = 5 # padding between frames
        py = 2 # padding between frames
        
        
        F1 = Frame(self.master)
        Label(F1, text='Geochemical Inputs').grid(row=0,columnspan=3)
        i = 1
        i = add_things_to_frame(F1, 10, 'Second Gas End-member', i)
        i = add_things_to_frame(F1, 11, 'Soil Gas End-member', i)
        i = add_things_to_frame(F1, 12, 'Mixed Gas', i)
        i = add_things_to_frame(F1, 17, 'Cave Air', i)
        F1.pack(side='left', fill=Y, anchor='n', padx=px, pady=py)
        
        F2 = Frame(self.master)
        Label(F2, text=' ').grid(row=0,columnspan=3)
        i = 1
        i = add_things_to_frame(F2, 16, 'Soil Metals (Chloride Salts)', i)
        i = add_things_to_frame(F2, 13, 'Bedrock Chemistry', i)
        i = add_things_to_frame(F2, 14, 'Bedrock Dissolution Conditions', i)
        i = add_things_to_frame(F2, 15, 'General', i)
        F2.pack(side='left', fill=Y, anchor='n', padx=px, pady=py)
        self.F2 = F2

        F3 = Frame(self.master)
        Label(F3, text='Model Scripting Options').grid(row=0,columnspan=3)
        i = 1
        i = add_things_to_frame(F3, 2, 'Scripting Options', i)
        i = add_things_to_frame(F3, 3, 'Additional PHREEQC output', i)
        i = add_things_to_frame(F3, 4, 'File IO Settings', i)

        # add run button
        RunButton = Button( F3, text="Run!", command= lambda : \
                            self._run_models())
        RunButton.grid(row=i, columnspan=3, sticky=S, pady=20)
        i = i + 1

        # add link to output GUI
        LinkButton = Button( F3, text="Open Output GUI", command= lambda : \
                             CCAnalyseGUI(Toplevel(self.master)) )
        LinkButton.grid(row=i, columnspan=3, sticky=S)
        
        F3.pack(side='left', anchor='n', padx=px, pady=py)
        self.F3 = F3
        
    def _run_models(self):
        
        s = self.settings.copy()
        out_dir = s.pop('out_dir').get()
        d = {}
        
        d1 = {k:v for (k,v) in s.items() if self.get_ln(k) != 'A'}
        d2 = {k:v for (k,v) in s.items() if self.get_ln(k) == 'A'}
        
        d1 = ns(tk2py(d1, parse=False))
        d2 = ns(tk2py(d2, parse=True))
        
        d = {**d1, **d2}

        p = cavecalc.forward_models.ForwardModels(settings=d, 
                                                  output_dir=out_dir)
        p.run_models()
        p.save()
        print("Done.")
        
class CCAnalyseGUI(object):
    """The Cavecalc Output GUI window."""
    
    def __init__(self, master):
        self.master = master
        self.master.title('Cavecalc Output GUI')
        self.e = Evaluate()
        self.dir = StringVar()
        self.dir.set(os.getcwd())
        self.loaded_dirs = []
        self.dnum = IntVar() # total no of models loaded
        self.dnum.set(0)
        self.settings_report = {}
        
        self.load_outputs_frame()
        self.save_buttons_frame()
        
    def _csv_dir_save(self):
        d = filedialog.askdirectory()
        if d:
            self.e.save_csvs(d)
        
    def _mat_save(self):
        d = filedialog.asksaveasfilename()
        if d:
            self.e.save_all_mat(file=d)
            
    def _add_data(self):
        """Loads data from the currently selected directory."""
        d = filedialog.askdirectory()
        
        if d:
            if d not in self.loaded_dirs:
                self.e.load_data(d)
                self.dnum.set(len(self.e.model_results))
                self.loaded_dirs.append(d)
        
    def load_outputs_frame(self):
        F0 = Frame(self.master)
        
        b = Button(master=F0, text="Load Model Output",
                   command = lambda : self._add_data())
        b.grid(row=0, column=0, columnspan=2)
        
        Label(F0, text="Models Loaded").grid(row=1, column=0, sticky=W)
        
        t1 = Entry(master=F0, textvariable=self.dnum, state='readonly')
        t1.grid(row=1,column=1)        
        F0.pack()
                  
    def save_buttons_frame(self):
        F1 = Frame(self.master)
        
        b1 = Button(master=F1, text="save as .csv", 
                    command = lambda : self._csv_dir_save())
        b1.grid(row=0, column=0)
        b2 = Button(master=F1, text="save as .mat", 
                    command = lambda : self._mat_save())
        b2.grid(row=0, column=1)
        b3 = Button(master=F1, text='Open Plotting Window',
                    command = lambda : PlottingWindow(self))
        b3.grid(row=0, column=2)
        
        F1.pack()
        
class PlottingWindow(Toplevel):
    def __init__(self, CCAnalyseGUI):
        super().__init__(CCAnalyseGUI.master)
        self.title('Cavecalc Plotting')
        
        self.e = CCAnalyseGUI.e
        self.o = self.e.model_results[0]
        self.s = self.e.model_settings[0]
        self.report = self.e.get_settings_report()
        
        lx = Label(self, text = "X variable (Required)")
        ly = Label(self, text = "Y variable (Required)")
        ll = Label(self, text = "Label with (Optional)")
        X_Sel, self.x = self.OutputSelectWidget()
        Y_Sel, self.y = self.OutputSelectWidget()
        L_Sel, self.l = self.SettingSelectWidget()
        
        lx.grid(row=0, column=0, sticky=W)
        X_Sel.grid(row=0, column=1)
        ly.grid(row=1, column=0, sticky=W)
        Y_Sel.grid(row=1, column=1)
        ll.grid(row=2, column=0, sticky=W)
        L_Sel.grid(row=2, column=1)
        
        # Radiobutton selector for data filtering (for plot)
        self.v = IntVar()
        self.v.set(0)
        r0 = Radiobutton(   self, text="Full Model (inc. initial solution)", 
                            variable=self.v, value=0)
        r1 = Radiobutton(   self, text="Full Model (excl. inital solution)", 
                            variable=self.v, value=1)
        r2 = Radiobutton(   self, text="Bedrock Dissolution Solution only", 
                            variable=self.v, value=2)
        r3 = Radiobutton(   self, text="End Point Solution only", 
                            variable=self.v, value=3)
        r4 = Radiobutton(   self, text="Precipitation Steps only",
                            variable=self.v, value=4)
                            
                            
        r0.grid(row=3)
        r1.grid(row=4)
        r2.grid(row=6)
        r3.grid(row=7)
        r4.grid(row=5)
                            
        b = self.PlotButton()
        b.grid(row=8,columnspan=2)

    def SettingSelectWidget(self):
        v = StringVar()
        opt = list(self.report.keys())
        if HIDDEN_OPTS:
            opt = [o for o in opt if o not in HIDDEN_OPTS]
        
        o = []
        for entry in opt:
            try:
                o.append(ns(entry))
            except KeyError:
                o.append(entry)
        return OptionMenu(self, v, *sorted(o)), v
        
    def OutputSelectWidget(self):
        v = StringVar()
        opt = list(self.o.keys())
        return OptionMenu(self, v, *sorted(opt)), v
    
    def PlotButton(self):

        b = Button(self, text='Plot Graph', 
                    command = lambda : self.plot())
        return b
        
    def plot(self):
    
        if self.v.get() == 1:
            a = self.e.filter_by_index(ind=0, n=True)
        elif self.v.get() == 2:
            a = self.e.filter_by_index(ind=1)
        elif self.v.get() == 3:
            a = self.e.filter_by_index(ind=-1)
        elif self.v.get() == 4:
            a = self.e.filter_by_results('step_desc', 'precip')
            
            
        else:
            a = copy.deepcopy(self.e)
           
        f = lambda x : None if x == '' else x
        x_lab = f(self.x.get())
        y_lab = f(self.y.get())
        lab_name = f(self.l.get())
        
        x = []
        y = []
        for v in a.model_results:
            x.append(v[x_lab])
            y.append(v[y_lab])
        
        if lab_name:
            labs = [s[ns(lab_name)] for s in a.model_settings]       
        else:
            labs = None
        
        print("Plotting...")
        gplot( x, y, x_lab, y_lab, labs, lab_name )
        
if __name__ == '__main__':
    CCInputGUI()
    # CCAnalyseGUI()
