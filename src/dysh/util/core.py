"""
Core utility classes and functions
"""

import numpy as np


def consecutive(data, stepsize=1):
    """Returns the indices of elements in `data`
    separated by less than stepsize separated into 
    groups.

    Parameters
    ----------
    data : array
        Array with values to split.
    stepsize : int
        Maximum separation between elements of `data`
        to be considered a single group.

    Returns
    -------
    groups : `~numpy.ndarray`
        Array with values of `data` separated into groups.
    """
    return np.split(data, np.where(np.diff(data) >= stepsize)[0]+1)


def sq_weighted_avg(a,axis=0,weights=None):
#@todo make a generic moment or use scipy.stats.moment
    r'''Compute the mean square weighted average of an array (2nd moment).

    :math:`v = \sqrt{\frac{\sum_i{w_i~a_i^{2}}}{\sum_i{w_i}}}`

    Parameters
    ----------
    a : `~numpy.ndarray`
        The data to average
    axis : int
        The axis over which to average the data.  Default: 0
    weights : `~numpy.ndarray` or None
        The weights to use in averaging.  The weights array must be the
        length of the axis over which the average is taken.  Default:
        `None` will use equal weights.

    Returns
    -------
    average : `~numpy.ndarray`
        The average along the input axis
    '''
    if weights is None:
        w = np.ones_like(a)
    else:
        w = weights
    v = np.sqrt(np.average(a*a,axis=axis,weights=weights))
    return v


def get_size(obj, seen=None):
    """Recursively finds size of objects.
    See https://goshippo.com/blog/measure-real-size-any-python-object/
    """
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size

def minimum_string_match(s,valid_strings):
    """return the valid string given a minimum string input"""
    pass

def stripTable(table):
    '''Remove leading and trailing chars from all strings from an input table.
        
       Parameters
       ----------
        table: ~astropy.table.Table
            The table to strip
    '''
    for n in table.colnames:
        if np.issubdtype(table.dtype[n],str):
            table[n] = np.char.strip(table[n])

#def strip(self,tables):
#    '''remove leading and trailing chars from all strings in list of tables'''
#    for b in tables:
#        stripTable(b)

def uniq(seq):
    """ Remove duplicates from a list while preserving order.
        from http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order 
    """
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

