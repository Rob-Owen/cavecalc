import numpy as np
import scipy.io as sio

def matlab_header_parse(dictionary):
    """Remove illegal characters from dictionary keys.
    
    To transfer data to matlab, certain characters are not allowed in field / 
    array names (e.g. brackets, hypens). This function removes these characters
    in preparation for writing dicts (of numpy arrays) to a .mat file.
    
    Args:
        dictionary: A dict
    Returns:
        A dict with modified key names
    """
    
    b = {}
    for k in dictionary:
        new_k = k.replace(')', '') # because matlab is fussy
        new_k = new_k.replace('(','_') # about variable names
        new_k = new_k.replace('-','') # and does not allow these
        new_k = new_k.replace('/','') # characters
        new_k = new_k.replace('[','')
        new_k = new_k.replace(']','')    
        
        b[new_k] = dictionary[k]
    return b
    
def numpify(dictionary):
    """Prepare a dict of lists for writing to a .mat file.
    
    Convert the dict of lists to a dict of numpy arrays. Dict keys are edited
    if they contain matlab-illegal characters. Lists in the dict are converted
    to numpy arrays.
    
    Args:
        dictionary: A dict of lists. Each list should be composed of a single
            type.
    Returns:
        A modified dict, ready for writing to a .mat file.    
    """
    
    a = matlab_header_parse(dictionary)
    b = {}
    for k in a:
        if k == "step_desc":
            b[k] = np.asarray( a[k], order='F' )
        else:
            b[k] = np.asarray( a[k], order='F' )    
    return b
            
def save_mat(dict_of_lists, filename):
    """Save data to a .mat file for use with Matlab/Octave.
    
    Takes a dict of lists (e.g. model results) and saves them to a .mat file.
    Data are prepared for saving using the numpify() function.
    
    Args:
        dict_of_lists: Data to be saved for Matlab use.
    """
        
    out = numpify(dict_of_lists)
    sio.savemat( filename, out )