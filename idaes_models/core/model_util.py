"""
This contains some general purpose functions that are commonly used in
various parts of IDAES models.
"""
from __future__ import print_function
from __future__ import division
__author__ = "John Eslick"
__version__ = "1.0.0"

def dict_set(v, d, pre_idx=None, post_idx=None, fix=False):
    """
    Set the values of array variables based on the values stored in a
    dictionary.  There may already be a better way to do this.  Should
    look into it.

    The value of Pyomo variable element with index key is set to d[key]

    Arguments:
    v: Indexed Pyomo variable
    d: dictonary to set the variable values from, keys should match a subset
        of Pyomo variable indexes.
    pre_idx: fixed indexes before elements to be set or None
    post_idx: fixed indexes after elements to be set or None
    fix: bool, fix the variables (otional)
    """
    #TODO: imporve doc string need to work out a good explaination <JCE>
    if pre_idx is None and post_index is None:
        for key in d:
            v[key].value = d[key]
            if fix:
                v[key].fixed = fix
    else:
        if pre_idx is None: pre_idx = ()
        if post_idx is None: post_idx = ()
        if not isinstance(pre_idx, tuple): pre_idx = (pre_idx,)
        if not isinstance(post_idx, tuple): post_idx = (post_idx,)
        for key in d:
            if not isinstance(key, tuple):
                key2 = (key,)
            else:
                key2 = key
            v[pre_idx + key2 + post_idx].value = d[key]
            if fix:
                v[pre_idx + key2 + post_idx].fixed = fix
