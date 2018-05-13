"""Generates settings for use in forward_models.py and caves.py for running
(suites of) models. SettingsObject instances are used for validating input
settings and passing them to caves.py to run.

Classes defined here:
    NameSwitcher
    SettingsObject
    SettingsMaker
"""

import os
import copy
import cavecalc.gui.mapping
import cavecalc.data.defaults as ccd
import cavecalc.data.types_and_limits as tal

cc_defaults = {k:v for k,v in vars(ccd).items() if '__' not in k}
cc_types = {k:v for k,v in vars(tal).items() if '__' not in k}
               
class NameSwitcher(object):
    """Handles switching parameter names between code and readable versions.
    
    Names and their readable equivalents are given in cavecalc.gui.mapping.
    
    Useage:
        ns = NameSwitcher()
        name2 = ns(name1)
        name3 = ns(name2)   # name3 = name1
    """
    
    def __init__(self):
        """Initialise a NameSwitcher object.
        
        Loads a name mapping ready for use.
        """
        
        self._name_map()
    
    def __call__(self, input):
        """Convert parameter names to their alternative formats.
        
        If names are passed as code names, readable versions are returned and
        vice versa.
        
        Args:
            input: Model name(s) to convert. May a string, list, set or dict
            where the keys represent the the parameter names.
        Returns:
            A copy of the original with the parameter names switched to their
            alternative forms.        
        """
        
        return self._switch_names(input)
    
    def _name_map(self):
        """Load name mappings as two dicts."""
        
        d1 = vars(cavecalc.gui.mapping)
        for k in copy.copy(d1):
            if '__' in k:
                d1.pop(k)
        d2 = {v : k for (k, v) in d1.items()}
        self.m2g = d1
        self.g2m = d2    
        
    def _switch_names(self, input):
        """Switches input entries (keys in the case of a dict) from model
        to gui parameter names, or vice versa.
        """
        
        s = input
        if type(s) is dict:
            try:
                output = { self.m2g[k] : v for (k,v) in s.items() }
            except KeyError:
                output = { self.g2m[k] : v for (k,v) in s.items() }
                    
        elif type(s) is set:
            try:
                output = { self.m2g[k] for k in s }
            except KeyError:
                output = { self.g2m[k] for k in s }
                
        elif type(s) is list:
            try:
                output = [ self.m2g[k] for k in s ]
            except KeyError:
                output = [ self.g2m[k] for k in s ]                    
        
        elif type(s) is str:
            try:
                output = copy.copy(self.m2g[s])
            except KeyError:
                output = copy.copy(self.g2m[s])
                
        else:
            raise TypeError("Input of type %s unsupported" % type(s))
        
        return output
        
class SettingsObject(object):
    """Constructs settings for use in Cavecalc models.
    
    Cavecalc caves.Simulator objects treats input settings as a dict and it is 
    possible to write these manually. However, use of SettingsObject is 
    advised - it provides default values, parameter name switching and value 
    checking.
    
    SettingsObjects are not normally instantiated directly - use SettingsMaker
    instead.
    """
    
    ns = NameSwitcher()
    
    def __init__(self, id=None, **kwargs):
        """Initialise settings object with default parameters.

        Default settings are stored in data/defaults.py.are loaded from 
        defaults.py. Non-default values may be specifed in kwargs, or later 
        using the set() method.
        
        Args:
            id: Model id number (optional)
            kwargs: Keyword arguments for model parameters
        """
        
        self.id = id
        
        # set defaults
        for k, v in cc_defaults.items():
            setattr(self, k, None)
            self.set(k, copy.copy(v))
            
        # set input parameters
        for k, v in kwargs.items():
            self.set(k, copy.copy(v))
        
    def __call__(self, *args):
        """Calls self.get() on *args."""
        
        return get(self, *args)
    
    def validate_entry(self, attribute, value='self'):
        """Check attribute, value pair is valid model input.
        
        Check entry against data/types_and_limits.py to confirm it is of 
        the correct type and/or value.
        
        Args:
            attribute (str): A model input parameter value.
            value: Value to be tested for validity. If not specified, the 
                value already set is used.
        """
        
        a = attribute
        
        if value == 'self':
            b = getattr(self, a)
        else:   
            b = value
        assert a in cc_types.keys()
        r = cc_types[a]

        ## for simple string or boolean input
        def type_test(t):
            def inner():
                if type(r) is t:
                    if isinstance(b, t):
                        return True
                    else:
                        raise TypeError("%s is of incorrect type: %r" %
                            (b, type(b)))
            return inner
        
        if type_test(str)(): return True
        if type_test(bool)(): return True
        
        ## for tuple-defined acceptable values
        if type(r) is tuple and isinstance(b, (float, int)):
            if r[0] is not None and b < r[0]:
                raise ValueError("%s invalid value: %f is less than %f" % 
                    (a, b, r[0]))
            if r[1] is not None and b > r[1]:
                raise ValueError("%s invalid value: %f is greater than %f" % 
                    (a, b, r[1]))
            return True         # numeric value within allowed range
        
        if type(r) is tuple and isinstance(b, str):
            try:
                if b in r[2:]:
                    return True # allowable string value detected
                else:
                    raise ValueError("%s invalid value: %r" % (a, b))
            except IndexError:
                raise TypeError("%s unsupported type: %r" % (a, type(b)))
                
        if type(r) is tuple:
            raise TypeError("%s unsupported type: %r" % (a, type(b)))
        
        ## for list-defined acceptable values
        if type(r) is list:
            if b.capitalize() in [e.capitalize() for e in r]:
                return True
            else:
                raise ValueError("%s invalid value: %r" % (a, b))
                
         
        raise Exception("Failed to check input type.\nData:\t%r\t%r\t%r"
                % (a, b, r))

    def set(self, parameter=None, value=None, **kwargs):
        """Set model input parameter(s).
        
        Arguments are first checked with validate_entry(). If validated, the
        model settings are updated.
        
        Args:
            parameter (str): The name of the parameter to be set.
            value: The value to set.
            **kwargs: Keyword arguments specifying multiple parameter=value
                pairs.
        """
        
        if parameter:
            kwargs[parameter] = value

        for k, v in kwargs.items():
            if k not in self.ns.m2g.keys():
                if k in self.ns.g2m.keys():
                    k = self.ns(k)
                else:
                    raise AttributeError(
                        "SettingsObject has no attribute %r" % k )
            
            if k == 'database':
                db_path = os.path.dirname(cavecalc.data.__file__)
                v = os.path.join(db_path,v)
            self.validate_entry(k, v)
            setattr(self, k, v)
    
    def get(self, *args):
        """Gets specified parameters from the object.
        
        If no arguments are specified, self.dict() is returned.
        
        Args:
            *args: names of parameters to be queried.
        Returns:
            A tuple of values corresponding to the length of *args. If args are
            not provided, self.dict() is returned.
        """
        
        if len(args) > 0:
            o = []
            for k in args:
                try:
                    o.append(getattr(self, k))
                except AttributeError:
                    o.append(getattr(self, self.ns(k)))
                    
            if len(o) == 1:
                return o[0]
            else:
                return tuple(o)
        else:
            return self.dict()
                
    def dict(self):
        """Return a dict of object settings.
        
        self.dict() is usually used to pass the model settings to a 
        caves.Simulator() object to be run. This excludes the 'id' setting, 
        which may be passed to Simulator separately.
        """
        
        a = copy.copy(vars(self))
        a.pop('id')
        return a
        
    def equals(self, SettingsObject):
        """Compare two SettingsObjects for equality.
        
        Args:
            SettingsObject: A second object to compare self to.
        Returns:
            True if all settings are equal. Otherwise False.
        """
        
        d1 = vars(self)
        d2 = vars(SettingsObject)
        
        mismatch = set(d1.items()) ^ set(d2.items())
        
        if len(mismatch) == 0:
            return True
        else:        
            return False

class SettingsMaker(object):
    """Generates suites of SettingsObjects.
    
    Allows easy generation of multiple SettingsObjects covering a range of 
    settings, usually to be run in series by cavecalc.forward_models."""
    
    def __init__(self, **kwargs):
        """Make a suite of SettingsObjects based on **kwargs.
        
        Args:
            **kwargs: keyword arguments of parameter=value for all non-default
            model input parameters. values may be single entries or lists of
            entries. If a lists are provided, a range of SettingsObjects are 
            generated, covering all possible combinations of input parameters.
            
        SettingsObjects generated are stored in the list self.o.
        """
        
        self.inps = kwargs
        self.sv = []    # single values (non-default constants)
        self.mv = []    # multiple values (a list of non-defaults)
        self.n = 1
        self.o = []

        self._sort_inps()
        
        if len(self.mv) > 0:
            self._make_set()
        else:
            self._make_one()
        self._set_ids()

        # self.print_report()
        
    def _sort_inps(self):
        """Split inputs into constants and variables"""
        
        for k,v in self.inps.items():
            if type(v) is list:
                self.mv.append((k, v))
            else:
                self.sv.append((k, v))

    def _make_one(self):
        consts_dict = {a[0] : a[1] for a in self.sv}
        self.o.append(SettingsObject( **consts_dict ))
        
    def _make_set(self, constants=None, variables=None):
        """Create a set of unique SettingsObjects covering all possible
        combinations of inputs. constants and variables are lists of key, 
        value tuples."""
    
        if not constants:
            constants = self.sv
        
        if not variables:
            variables = self.mv
                    
        if len(variables) > 1:
            a_var = variables.pop()
            for e in a_var[1]:
                consts = copy.deepcopy(constants)
                consts.append((a_var[0], e))
                self._make_set( consts, copy.deepcopy(variables) )
            return
            
        if len(variables) == 1:
            for e in variables[0][1]:
                consts_dict = {a[0] : a[1] for a in constants}
                variables_dict = {variables[0][0] : e}
                self.o.append(SettingsObject( **consts_dict, **variables_dict ))
                
    def _remove_copies(self):
        """Filters self.o to remove identical SettingsObjects."""
        
        def all_diff(lst):
            for i,a in enumerate(lst):
                for j,b in enumerate(lst):
                    if a.equals(b) and i!=j:
                        return False
            return True
    
        def del_ifn_unique():
            c = copy.copy(self.o)
            for i,a in enumerate(c):
                for j,b in enumerate(c):
                    if a.equals(b) and i!=j:
                        self.o.pop(i)
                        return
            
        while not all_diff(self.o):
            for i,m in enumerate(self.o):
                del_ifn_unique()
        
    def _set_ids(self):
        for i,a in enumerate(self.o):
            a.id = i
    
    def print_report(self):
        """Debug function for viewing SettingsObjects in self.o."""
        
        def list_cmp(l):
            for i,e in enumerate(l):
                try:
                    if e != l[i+1]:
                        return True
                except:
                    return False
            
        d = self.o[0].dict()
        for k in d.keys():
            v = []
            for m in self.o:
                v.append(m.get(k))
            if list_cmp(v):
                print(k,end='\t\t')
                for e in v:
                    print(e,end='\t')
                print('')
 
    def settings(self):
        """Return the list of generated SettingsObjects."""
        
        return copy.deepcopy(self.o)
 
if __name__ == '__main__': # testing code
    # a = SettingsObject()
    # b = SettingsObject()
    # print("Equal?" + '\t' + str(a.equals(a)))
    # print("Equal?" + '\t' + str(a.equals(b)))
    # a.set(temperature=35)
    # print("Equal?" + '\t' + str(a.equals(b)))
    # print(a.get('temperature'))
    
    # n = NameSwitcher()
    # print(n('temperature'))
    
    a = SettingsMaker(temperature=[15, 20], soil_pCO2=[2222, 3333])
    a.print_report()
    a = SettingsMaker()
    
    