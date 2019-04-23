"""
Functions for saving and loading Pyomo objects to json
"""
# Changes the divide behavior to not do integer division
from __future__ import division
from __future__ import print_function

from pyomo.environ import *
from pyomo.dae import *
from collections import OrderedDict
import six
import json
from six import itervalues, iteritems

# Some more inforation about this module
__author__ = "John Eslick <john.eslick@netl.doe.gov>"
__version__ = "1.0.0"


def _write_dict(o, sd=None):
    """
    Write supported Pyomo object state to a dictionary which can be exported
    to json.

    Args:
    o: Pyomo object to save
    sd: dict to add on to, if none just return new dict
    """
    if isinstance(o, Block) and not isinstance(o, pyomo.core.base.block._BlockData):
        d = _write_indexed_block_dict(o)
    elif isinstance(o, ConcreteModel):
        d = _write_model_dict(o)
    elif isinstance(o, Block):
        d = _write_block_dict(o)
    elif isinstance(o, pyomo.core.base.var.IndexedVar):
        d = _write_indexed_var_dict(o)
    elif isinstance(o, DerivativeVar):
        d = _write_derivative_var_dict(o)
    elif isinstance(o, Var):
        d = _write_var_dict(o)
    elif isinstance(o, pyomo.core.base.constraint.IndexedConstraint):
        d = _write_indexed_constraint_dict(o)
    elif isinstance(o, Constraint):
        d = _write_constraint_dict(o)
    else: # An object type that we don't currently save.
        return None
    # Return dict
    if sd is not None: # and not adding to existing dictionary
        sd[o.name] = d
    return d

def _read_dict(o, d, **kwargs):
    """
    Read a dict and set properties of o.
    """
    if d is None:
        pass # d may be none if object in the hierachy that's not supported
    elif isinstance(o, Block) and not isinstance(o, pyomo.core.base.block._BlockData):
        _read_indexed_block_dict(o, d, **kwargs)
    elif isinstance(o, ConcreteModel):
        _read_model_dict(o, d, **kwargs)
    elif isinstance(o, Block):
        _read_block_dict(o, d, **kwargs)
    elif isinstance(o, DerivativeVar):
        _read_derivative_var_dict(o, d, **kwargs)
    elif isinstance(o, pyomo.core.base.var.IndexedVar):
        _read_indexed_var_dict(o, d, **kwargs)
    elif isinstance(o, Var):
        _read_var_dict(o, d, **kwargs)
    elif isinstance(o, pyomo.core.base.constraint.IndexedConstraint):
        pass
    elif isinstance(o, Constraint):
        pass

def _write_model_dict(blk):
    """
    Write a dict to store data in an Pyomo concerete model block
    """
    b = {}
    d = {'metadata': {'type': 'model'}, 'data': b, 'active':blk.active}
    for o in blk.component_objects(descend_into=False):
        _write_dict(o, b)
    return d

def _read_model_dict(blk, sd, **kwargs):
    """
    Load a Pyomo ConcreteModel
    """
    sd = sd['data']
    for o in blk.component_objects(descend_into=False):
        d = sd.get(o.name, None)
        # if d is none is an object that wasn't saved (is okay)
        _read_dict(o, d, **kwargs)

def _write_indexed_block_dict(blk):
    """
    Write a dict to store data in an Pyomo indexd block
    """
    ib = OrderedDict()
    d = {'metadata': {'type': 'indexed_block'}, 'data': ib, 'active':blk.active}
    for i, b in blk.iteritems():
        ib[i] = _write_block_dict(b)
    return d

def _read_indexed_block_dict(blk, sd, **kwargs):
    """
    Load a Pyomo indexed block
    """
    sd = sd['data']
    dl = sd.values()
    if six.PY3:
        dl = list(dl)
    for i, ii in enumerate(blk):
        _read_dict(blk[ii], dl[i], **kwargs)

def _write_block_dict(blk):
    """
    Write a dict to store data in an Pyomo single block
    """
    b = {}
    d = {'metadata': {'type': 'single_block'}, 'data': b, 'active':blk.active}
    for o in blk.component_objects(descend_into=False):
        _write_dict(o, b)
    return d

def _read_block_dict(blk, sd, **kwargs):
    """
    Load a Pyomo single block
    """
    sd = sd['data']
    for o in blk.component_objects(descend_into=False):
        d = sd.get(o.name, None)
        # if d is none is an object that wasn't saved (is okay)
        _read_dict(o, d, **kwargs)

def _write_indexed_var_dict(var):
    """
    Write a dict to store data in an Pyomo indexed variable
    """
    iv = OrderedDict()
    d = {'metadata': {'type': 'indexed_var'}, 'data': iv, 'active':var.active}
    for i, v in var.iteritems():
        iv[str(i)] = _write_var_dict(v)
    return d

def _read_indexed_var_dict(var, sd, **kwargs):
    """
    Load Pyomo indexed variable
    """
    sd = sd['data']
    dl = sd.values()
    if six.PY3:
        dl = list(dl)
    for i, ii in enumerate(var):
        _read_dict(var[ii], dl[i], **kwargs)

def _write_derivative_var_dict(var):
    """
    Write a dict to store data in an Pyomo derivative variable
    """
    iv = OrderedDict()
    d = {'metadata': {'type': 'derivative_var'}, 'data': iv, 'active':var.active}
    for i, v in var.iteritems():
        iv[str(i)] = _write_var_dict(v)
    return d

def _read_derivative_var_dict(var, sd, **kwargs):
    """
    Load Pyomo derivative var
    """
    sd = sd['data']
    dl = sd.values()
    if six.PY3:
        dl = list(dl)
    for i, ii in enumerate(var):
        _read_dict(var[ii], dl[i], **kwargs)

def _write_var_dict(var):
    """
    Write a dict to store data in an Pyomo single variable
    """
    v = {}
    d = {'metadata': {'type': 'single_var'}, 'data': v, 'active':var.active}
    v['ub'] = var.ub
    v['lb'] = var.lb
    v['value'] = var.value
    v['fixed'] = var.fixed
    return d

def _read_var_dict(var, sd, **kwargs):
    """
    Load Pyomo single variable
    """
    value_only = kwargs["value_only"]
    bounds_only = kwargs["bounds_only"]
    sd = sd['data']
    try:
        if not bounds_only:
            var.value = sd['value']
    except:
        print("No value: {0}".format(var.cname()))
        print(sd)
    if not value_only:
        var.setub(sd['ub'])
        var.setlb(sd['lb'])
        if not bounds_only:
            if sd['fixed']:
                var.fix()
            else:
                var.unfix()

def _write_indexed_constraint_dict(con):
    ic = OrderedDict()
    d = {'metadata':{'type': 'indexed_constraint'},
         'data': ic,
         'active':con.active}
    for i, c in con.iteritems():
        ic[str(i)] = _write_constraint_dict(c)
    return d

def _write_constraint_dict(con):
    d = {'metadata':{'type': 'constraint'},
         'data': {},
         'active':con.active}
    return d

def save_json(o, fname=None, dict_out=False, human_read=False):
    """
    Saves Pyomo variable states. The result of this is a json or dict
    representation of a Pyomo object and hierarchy.  Each object consists of 2
    enrties metadata and data.  Metadata contains information about the object
    like what type it is.  Data contains all the Pyomo stuff being saved.

    Ordered dicts are used to keep indexed components in order without
    needing to worry too much about what the indexes acually are. The
    indexes are also saved as strings, but they may be harder to parse.
    Tuple indexes are not readily json serializable, but relying on
    maintaining the order of the elements should work well. So don't
    forget to use "object_pairs_hook=OrderedDict" when reading the json
    back or everything will be mixed up.

    Single block and single var are elements in indexed vars or blocks

    Args:
    o: A Pyomo object to save (also saves sub-objects)
    fname: A file to save the json as, if None don't save files
    dict_out: True, return Python dict; False, return json string
              dict output is useful for copying values between like
              blocks or other model substructures or models
    human_read: is mostly for debuging, makes the json easier to read
    """
    d = _write_dict(o)
    if dict_out:
        # if dict_out just return a Pyomo dict.  This is usually for
        # copying values between like objects
        return d
    # Else carry on and write json
    dump_kw = {'indent': 2} if human_read else {'separators': (',', ':')}
    if fname:
        # with file, write JSON and return None
        with open(fname, 'w') as ofile:
            json.dump(d, ofile, **dump_kw)
        return
    else:
        # with no file, return JSON string
        return json.dumps(d, **dump_kw)

def load_json(o, fname=None, s=None, sd=None, **kwargs):
    """
    Loads Pyomo object variable values from JSON.  See the
    corresponding save_json function.

    Args:
    fname: json file path to load, if None load from differnt thing
    s: A json string to load, if None load from a differnt thing
    sd: A Python dict to load, if None load a differnt object
    values_only: Just load variable values and nothing else.

    Must specify exactly one of fname, s, or sd.
    """
    kwargs.setdefault("value_only", False)
    kwargs.setdefault("bounds_only", False)
    # Should only one but not checking
    if sd is not None:
        d = sd
    elif fname is not None:
        with open(fname, 'r') as f:
            d = json.load(f, object_pairs_hook=OrderedDict)
    elif s is not None:
        d = json.loads(s, object_pairs_hook=OrderedDict)
    else:
        raise Exception("Need to specify one source.")
    # Read data
    _read_dict(o, d, **kwargs)
