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

#Special callback unctions for reading or writing attibutes
#that cannot be set or read directly.
def _set_active(o, d):
    """
    Set if component is active
    """
    if d:
        o.activate()
    else:
        o.deactivate()

def _set_fixed(o, d):
    """
    Set if variable is fixed
    """
    if d:
        o.fix()
    else:
        o.unfix()

def _set_lb(o, d):
    """
    Set variable lower bound
    """
    o.setlb(d)

def _set_ub(o, d):
    """
    Set variable upper bound
    """
    o.setub(d)


class StoreSpec(object):
    """
    Object to tell read/write json what to read or write and how
    """
    def __init__(self,
                 classes=(
                    (Param, ()),
                    (Var, ()),
                    (Component, ("active",))),
                data_classes=(
                    (pyomo.core.base.var._VarData,
                        ("fixed", "stale", "value", "lb", "ub")),
                    (pyomo.core.base.param._ParamData,
                        ("value",)),
                    (pyomo.core.base.component.ComponentData, ("active",))),
                skip_classes=
                    (ExternalFunction, Set, Connector, Expression, RangeSet),
                tree_struct=True):
        """
        Initialize an object to specify what parts of a model are saved.
        classes and data classes are checked in order.  So the more specific
        classes should go first and fallback cases should go last.  Since
        classes like component catch pretty muh everything, you can also
        specify classes to skip.

        Args:
        classes: The classes to save, and what attributes to save, if the
            attributes are not directly readable/writable you can set callbacks
            for them later.
        data_classes: The component data classes to save, and what attributes to
            save, if the attributes are not directly readable/writable you can
            set callbacks for them later.
        skip_classes: Component classes to explicity skip reading or writing
        tree_struct: save the parent and child info.  the stucture of this file
            doesn't preserve the hierarchy otherwise (other that could be
            deduced from the component names).
        """
        self.write_cbs={}
        self.read_cbs={
            "active":_set_active,
            "fixed":_set_fixed,
            "lb":_set_lb,
            "ub":_set_ub}
        self.skip_classes = skip_classes
        self.classes = [i[0] for i in classes]
        self.data_classes = [i[0] for i in data_classes]
        self.class_attrs = [i[1] for i in classes]
        self.data_class_attrs = [i[1] for i in data_classes]
        self.tree_struct = tree_struct

    def set_read_callback(self, attr, cb=None):
        """
        Set a callback to set an attribute, when reading from json or dict.
        """
        self.read_cbs[attr] = cb

    def set_write_callback(self, attr, cb=None):
        """
        Set a callback to get an attribute, when writing to json or dict.
        """
        self.write_cbs[attr] = cb

def StoreSpecBounds():
    return StoreSpec(
        classes=((Var, ()),),
        data_classes=((pyomo.core.base.var._VarData, ("lb", "ub")),),
        tree_struct=False)

def _write_component(sd, o, wts):
    """
    Writes a componet to the save dict
    """
    for i, cl in enumerate(wts.skip_classes):
        if isinstance(o, cl):
            return
    alist = None
    for i, cl in enumerate(wts.classes):
        if isinstance(o, cl):
            alist = wts.class_attrs[i]
            break
    if alist is None: return
    sd[o.name] = {}
    if wts.tree_struct:
        if o._parent is None:
            sd[o.name]["parent"] = None
        else:
            parent=o._parent()
            sd[o.name]["parent"] = parent.name
    for a in alist:
        d = getattr(o, a, None)
        sd[o.name][a] = d
    #Write data elements
    alist = []
    sd[o.name]["data"] = {}
    indexed=True
    c=0
    for key, el in o.iteritems():
        if c == 0:
            for i, cl in enumerate(wts.data_classes):
                if isinstance(el, cl):
                    alist = wts.data_class_attrs[i]
                    break
        c+=1
        if key is None or key=="": indexed=False
        sd[o.name]["data"][str(key)] = {}
        for a in alist:
            d = getattr(el, a)
            sd[o.name]["data"][str(key)][a] = d
    if c>1: indexed=True
    sd[o.name]["indexed"] = indexed

def save_json(o, fname=None, human_read=False, wts=None):
    if wts is None:
        wts = StoreSpec()
    sd={}
    _write_component(sd, o, wts)
    for o in o.component_objects(descend_into=True):
        _write_component(sd, o, wts)
    if fname is not None:
        dump_kw = {'indent': 2} if human_read else {'separators': (',', ':')}
        with open(fname, "wb") as f:
            json.dump(sd, f, **dump_kw)
    return sd

def _read_component(sd, o, wts):
    """
    Read a component dictionary into a model
    """
    for i, cl in enumerate(wts.skip_classes):
        if isinstance(o, cl):
            return
    alist = None
    for i, cl in enumerate(wts.classes):
        if isinstance(o, cl):
            alist = wts.class_attrs[i]
            break
    if alist is None: return
    for a in alist:
        if a in wts.read_cbs:
            wts.read_cbs[a](o, sd[o.name][a])
        else:
            setattr(o, a, sd[o.name][a])
    #Write data elements
    alist = []
    c = 0
    for key, el in o.iteritems():
        if c==0:
            for i, cl in enumerate(wts.data_classes):
                if isinstance(el, cl):
                    alist = wts.data_class_attrs[i]
                    break
        c +=1
        for a in alist:
            if a in wts.read_cbs:
                wts.read_cbs[a](el, sd[o.name]["data"][str(key)][a])
            else:
                setattr(el, a, sd[o.name]["data"][str(key)][a])

def load_json(o, sd=None, fname=None, s=None, wts=None):
    """
    load a json file of python dict
    """
    if sd is not None:
        pass #python dict
    elif fname is not None:
        with open(fname, "rb") as f:
            sd=json.load(f) #json file
    elif s is not None:
        sd=json.loads(s) #json string
    else:
        raise Exception("Need to specify a data source to load from")
    if wts is None:
        wts = StoreSpec()
    _read_component(sd, o, wts)
    for o in o.component_objects(descend_into=True):
        _read_component(sd, o, wts)
