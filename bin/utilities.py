import time 
import logging

def timing(f):
    """
    Helper function for timing other functions
    Parameters
    ----------
    f : function
    Returns
    -------
    function
        new function wrap with timer and logging 
    """
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        logging.debug('{:s} function took {:.10f} s'.format(f.__name__, (time2-time1)))
        return ret
    
    return wrap