"""
Base for IDAES process model objects.
"""
from __future__ import division  # No integer division
from __future__ import print_function  # Python 3 style print

from pyomo.environ import *
from pyomo.dae import *
from pyomo.core.base.block import SimpleBlock, IndexedBlock, _BlockData
from pyomo.core.base.external import AMPLExternalFunction
from pyomo.core.base.connector import ConnectorExpander
# import copy
import os
import sys
from datetime import datetime
import networkx as nx
import weakref
import itertools
import pyomo.environ as pe
from six import itervalues, iteritems, string_types
from bunch import Bunch
from .util.oa import apply_OA, add_oa_constraints
from .util.solve import SolverError, get_solvers
from .util.idjson import save_json, load_json
from .util.var import Wrapper
#from pyomo.core.base.indexed_component import _IndexedComponent_slice
from .util import idjson2 as j2
# Some more inforation about this module
__author__ = "John Eslick <john.eslick@netl.doe.gov>, "\
             "Qi Chen <qichen@andrew.cmu.edu>"
__version__ = "1.0.0"

try:
    import matplotlib.pyplot as plt
except:
    plt = None

__all__ = ['ProcessBase', 'ProcessBlock']

# Reserved keywords that go to the Block constructor 'name' is also a block
# keyword, but I'm stealing it, so it's not on the list
_block_kwds = ('rule', 'options', 'concrete', 'ctype', 'noruleinit', 'doc')

def ProcBlock(name):
    def ProcBlockDec(cls):
        c = type(name, (ProcessBlock,), {"_data_class":cls, "__module__":cls.__module__})
        setattr(sys.modules[cls.__module__], name, c)
        return cls
    return ProcBlockDec


class _InitMethod():
    """Class contains a function and documentation string for an intialization
    method.
    """

    def __init__(self, func, doc=""):
        """
        Args:
        func: a function to perform a model intialization routine
        doc: a description of the initialization routine
        """
        self.func = func
        self.doc = doc
        self.long_doc = ""


class CallableDict(dict):
    """
    A dictionary that returns itself when called.  Behaves like a weakref would.

    The reason I have this is so I can treat a dictionary containing Pyomo vars
    the same as a weakref to an indexed Pyomo var.  It helps with dealing with
    links to Pyomo variable slices
    """
    def __call__(self):
        return self


class ExternalFunction(AMPLExternalFunction):
    """
    Workaround for a bug in Pyomo.  DAE transformation complains of missing
    dim() method for ExternalFunction objects the -1 return value makes it
    dodge some stuff in the DAE transoframtion code.  I'm assuming dim isn't
    going to cause problems elsewhere.  DAE was the only thing that missed it.
    """
    # TODO: remove when Pyomo is fixed. <JCE>
    def dim(self):
        return -1

def _pop_nonblock(dct):
    """
    Split a kwargs dict into kwargs meant for block and args meant for the thing
    that inherits block data.
    """
    da = {}
    pop_list = [key for key in dct if key not in _block_kwds]
    for key in pop_list:
        da[key] = dct.pop(key, {})
    return da


class _IndexedProcessBlockMeta(type):
    """
    This is a metaclass used to create an indexed model class
    """
    def __new__(meta, name, bases, dct):
        def __init__(self, *args, **kwargs):
            kwargs.pop("data_class", None)
            _pop_nonblock(kwargs)
            Block.__init__(self, *args, **kwargs)
        dct["__init__"] = __init__
        return type.__new__(meta, name, bases, dct)


class _SimpleProcessBlockMeta(type):
    """
    This is a metaclass used to create a simple model class
    """
    def __new__(meta, name, bases, dct):
        def __init__(self, *args, **kwargs):
            kwargs.pop("data_class", None)
            da = _pop_nonblock(kwargs)
            bases[0].__init__(self, self, idx=None, **da)
            kwargs["concrete"] = True
            Block.__init__(self, *args, **kwargs)
            self.__init2__()
            self._data[None] = self
        dct["__init__"] = __init__
        return type.__new__(meta, name, bases, dct)


class ProcessBlock(Block):
    """
    """
    def __new__(cls, *args, **kwds):
        if cls.__name__.startswith('_Indexed') or \
            cls.__name__.startswith('_Simple'):
            return super(Block, cls).__new__(cls)
        _data_args = _pop_nonblock(kwds)
        if args == ():
            bname = "_Simple{}".format(cls.__name__)
            n = _SimpleProcessBlockMeta(bname, (cls._data_class, cls), {})
            return n.__new__(n)
        else:
            bname = "_Indexed{}".format(cls.__name__)
            n = _IndexedProcessBlockMeta(bname, (cls,), {})
            o = n.__new__(n)
            o._data_class = cls._data_class
            o._data_args = _data_args
            return o

    def _default(self, idx):
        return self._data.setdefault(
            idx, self._data_class(self, idx=idx, **self._data_args))

    def save_json(self, *args, **kwargs):
        """
        Save Pyomo object state to JSON.  See .util.idjson import save_json
        and load_json for more details.
        """
        return j2.save_json(self, *args, **kwargs)

    def load_json(self, *args, **kwargs):
        """
        Loads Pyomo object state from json.  See .util.idjson import save_json
        and load_json for more details.
        """
        j2.load_json(self, *args, **kwargs)


class ProcessBase(_BlockData):
    """
    This is the IDAES model base class for process model objects like
    unit models, property models, reaction models, ...

    These model objects construct and manage Pyomo Blocks or Indexed
    Blocks

    Attributes:
        bstore (TYPE): Description
        con_graph (TYPE): Description
        constraint_types (dict): Dictionary of constraint types. Contains two
            elements: 'linear' and 'nonlinear'.
        delay_construct (Boolean): flag for whether construction of the class
            should be delayed.
        fluid_connector_elements (dict): Description
        graph (TYPE): Description
        link_map (dict): Dictionary storing link information between sub-units
            contained within this class.
        links (Block): Pyomo Block containing linking constraints between
            sub-units.
        model_type (str): Description
        process_graph (TYPE): Description
        pyomo_results (TYPE): Description
        pyomo_solver (TYPE): Description
        solve_time (TYPE): Description
        solvers (dict): dictionary of solvers for various subproblems
        units (dict): Dictionary of process units contained in this class.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize the object.  Anything inheriting from process base
        should impliment a build() function that create the Pyomo
        variables, constriants, and whatever.

        Args:
        ------
        delay_construct: A model that inherits this can call the build
            function at the end of the constructor (or not).  The
            delay_construct argument provides a standard way to override
            the automatic call to build and manually do it later. The
            purpose of delaying the call to build is to do some more
            model setup before constructing the Pyomo problem. False by
            default
        solver: the solver to be used in solving the this Pyomo Block
            by default is ipopt.  The solver and solver settings can
            be changed at anytime.
        """
        # Pop arguments that are used here and not to be passed on
        # to the Pyomo Block constructor.
        kwargs.setdefault('name', 'Unnamed_Unit_Model')
        self.delay_construct = kwargs.pop("delay_construct", False)
        self.unit_type = kwargs.pop('type', ['General'])
        self.unit_name = kwargs.pop('name')
        self.parent_unit = Wrapper(kwargs.pop('parent'))

        self._idx = kwargs.pop('idx', None)
        self._comps = set(kwargs.pop('comps', []))

        self._init_methods = {}
        self.link_map = {}
        self.units = Bunch()
        self.graph = nx.DiGraph()
        self._built = False
        _BlockData.__init__(self, *args)
        #super(ProcessBase, self).__init__(self, *args)
        if self.parent_component().__class__.__name__.startswith("_Indexed"):
            #if its a simple block this gets called from Simple* init.  need
            #Block init in between __init__ and __init2__
            self._suppress_ctypes = self.parent_component()._suppress_ctypes
            self.__init2__()
        self._reg_initialize(
            'default',
            func=self.default_init,
            doc='Do nothing')

    def __init2__(self):
        """
        Second part of init to add Pyomo objects, which can't be done until after
        block init
        """
        if isinstance(self.unit_type, string_types):
            self.unit_type = [self.unit_type]
        if self.parent_unit._obj is None:
            self.parent_unit = None  # do not wrap None
        self.solvers = get_solvers()
        self.links = pe.Block()
        self.fetch_limits()
        if not self._comps:
            try:
                self._comps = self.parent_unit._comps
            except:
                pass
        self.comps = pe.Set(initialize=self._comps)
        if not self.delay_construct:
            self.build()

    def add_weakref_links(self, **kwargs):
        """
        This adds a weakref object pointing to a variable contained in another
        block.
        """
        for key in kwargs:
            #expand slices
            o = kwargs[key]
            if isinstance(o, _IndexedComponent_slicer):
                indexes = []
                indexes2 = []
                objs = []
                for o2 in o:
                    indexes.append(o2.index())
                    objs.append(o2)
                mask = [False]*len(indexes[0])
                for j in range(len(indexes[0])):
                    for i in indexes:
                        if indexes[0][j] != i[j]:
                            mask[j] = True
                            break
                d = CallableDict()
                setattr(self, key, d)
                for j, i in enumerate(indexes):
                    i = tuple(itertools.compress(i, mask))
                    if len(i) == 1:
                        i = i[0]
                    d[i] = objs[j]
            else:
                setattr(self, key, weakref.ref(o))

    def var_or_expr(self, *args, **kwargs):
        """
        Creates either an expression or an expression, variable, and constraints
        this allows the problem to be either constructed with more variables
        and simpler constraints or less variables and more complex constraints.
        If using expression only, and you want the find the values of some
        quantity after solving a model the expression can be evaluated.  If using
        vaiables the expression can be evaluated in an initialization proceedure.

        The arguments are the same as Var except:

        Args
        name: Name of variable or exprssion to create, if variable an
            Expression object (expr_name) and a constraint object (eq_name)
            will also be created.
        expr or rule must also be specified they are the same as the Expression
            expr and rule arguments
        """
        if "domain" not in kwargs:
            kwargs["domain"] = Reals
        name = kwargs.pop("name")
        expr_only = kwargs.pop("expr_only", self.expr_only)
        eq_scale = kwargs.pop("eq_scale", 1.0)
        if "rule" in kwargs:
            expr = Expression(*args, rule=kwargs.pop("rule"))
        elif "expr" in kwargs:
            expr = Expression(*args, expr=kwargs.pop("expr"))
        else:
            raise Exception("expr or rule required")
        if expr_only:
            setattr(self, name, expr)
        else:
            setattr(self, "expr_"+name, expr)
            if "initialize" not in kwargs and not args:
                kwargs["initialize"] = value(expr)
            v = Var(*args, **kwargs)
            setattr(self, name, v)
            if "initialize" not in kwargs and args:
                for i in v:
                    v[i].value = value(expr[i])
            if not args:
                setattr(self, "eq_"+name,
                    Constraint(expr=eq_scale*v==expr*eq_scale))
            else:
                def eqrule(blk, *index):
                    return eq_scale*v[tuple(index)]==eq_scale*expr[tuple(index)]
                setattr(self, "eq_"+name, Constraint(*args, rule=eqrule))

    def default_init(self):
        """
        An model initialization function that does nothing.
        """
        pass

    def unfix_all(self):
        """This unfixes all variables in this block and sub blocks. This is
        mostly for initialization to help get the model into a known state to
        start the initilaization proceedure.
        """
        pass

    def activate_all(self):
        """This acivates all Pyomo objcects in this block and subblocks. This is
        mostly for initialization to help get the model into a known state to
        start the initilaization proceedure.
        """
        pass

    def load_external_funcs(self, lib_file, func_dict):
        """
        Load external functions that can be used in Pyomo expressions

        Args:
        lib_file: A shared library file (so, dll...) containing ASL user defined
            functions
        func_dict: Dictionary of functions to load, key = attribute name,
            value = name of user ASL function in lib_file
        """
        for f in func_dict:
            setattr(self, f, ExternalFunction(library=lib_file,
                    function=func_dict[f]))

    def _fluid_connector_constraints(self):
        """This returns a dictionary of constraints that connect two
        fluid connectors and a list of constraints that result from
        their expansion.
        """
        condict = {}
        explist = []
        for o in self.component_objects(descend_into=True):
            if isinstance(o, Constraint):  # find all constraints
                if o.name[-9:] == ".expanded":
                    explist.append(o)
                if getattr(o, "connection", False) == "fluid":
                    condict[o.name] = o
        return (condict, explist)

    def fluid_connection_deactivate(self, name):
        """This deactivates constraints that result from the expation of
        a constraint (named name) equating to fluid connectors.
        """
        self.fluid_connection_activate(name=name, active=False)

    def push_connection(self, name):
        """
        Set the values in a 'to' connector to the values in a from connector.
        This helps with sequential modular initialization by move values from
        outlets of one unit to inlet
        """
        lnk = self.link_map[name]
        c1 = getattr(lnk["from"], lnk["from_port"])
        c2 = getattr(lnk["to"], lnk["to_port"])
        for vkey in c1.vars:
            v = c1.vars[vkey]
            if isinstance(v, pyomo.core.base.var.IndexedVar):
                for i, ve in v.iteritems():
                    try:
                        c2.vars[vkey][i].value = ve.value
                    except:
                        raise("Pushing vales along incompatible connetor")
            else:
                try:
                    c2.vars[vkey].value = v.value
                except:
                    raise("Pushing vales along incompatible connetor")

    def fix_connector(self, name, fixed=True):
        """
        Fix or unfix the variables in a connector.  This is used mainly to make
        sequential modular type initilialization a little cleaner.

        Args:
        name: name of the connector
        fixed: Frue to fix, False to unfix
        """
        c = getattr(self, name, None)
        if c is None:
            raise Exception("Connector '{}' does not exist".format(name))
        elif not isinstance(c, Connector):
            raise Exception("'{}' is not a connector".format(name))

        for vkey in c.vars:
            v = c.vars[vkey]
            if v is None: continue
            if isinstance(v, pyomo.core.base.var.IndexedVar):
                for i, ve in v.iteritems():
                    if fixed:
                        ve.fix()
                    else:
                        ve.unfix()
            else:
                if fixed:
                    v.fix()
                else:
                    v.unfix()

    def fluid_connection_activate(self, name, active=True):
        """This acivates or deactivates a constraints that result from
        the expation of a constraint (named name) equating to fluid
        connectors.
        """
        condict, explist = self._fluid_connector_constraints()
        for o in explist:
            if o.name[:-9] == name:
                if active:
                    o.activate()
                else:
                    o.deactivate()

    def add_fluid_port(self, *args, **kwargs):
        """ Currently adds a fluid connector in the future it may also
        add a port or a connector based on the multiport flag.

        Reserving args for creating indexed connectors, currently not
        used, but may be a future addition. All arguments must be kwargs

        args:
        single_choice: exclusivly one inlet or outlet for multiport
        multi_port: allow multiple possible inlet/outlet connections
        direction: indicate the port's intended flow direction
        F: flow variable from unit model to be connected
        fc: component flow variable
        T: temperature variable
        P: pressure variables
        vf: vapor fraction variable
        y: mole fraction variable assume I can get component list here
        """
        name = kwargs.pop("name", None)
        if name is None:
            raise Exception("name argument required")
        # single_choice = kwargs.pop("single_choice", False)
        # multi_port = kwargs.pop("multi_port", False)
        direction = kwargs.pop("direction", "in")
        F = kwargs.pop("F", None)
        fc = kwargs.pop("fc", None)
        T = kwargs.pop("T", None)
        P = kwargs.pop("P", None)
        vf = kwargs.pop("vf", None)
        y = kwargs.pop("y", None)
        # add connector
        c = Connector()
        c.direction = direction
        c.port_type = "fluid"
        # c.direction = direction
        c.add(F, name="F")
        c.add(fc, name="fc")
        c.add(T, name="T")
        c.add(P, name="P")
        c.add(y, name="y")
        c.add(vf, name="vf")
        setattr(self, name, c)

        # Leave this to rember, when building a graph can just look into
        # the components, see if they are connectors and
        # maybe get partent if needed.  Dont need to keep track in a dict
        # just reconstruct as needed?
        # print(c._parent)
        # print(c.name)

    def add_destination(self, to_port, from_port=''):
        """Registers a destination unit port 'to_port' with the outlet port of
        this unit denoted by 'from_port'

        Args:
            to_port (UnitPort): the port on the destination unit
            from_port (str, optional): the name of the port on this unit to
                associate with the destination unit port. If not specified,
                then a default name of 'outlet' is used.

        Returns:
            None
        """
        from_port = 'outlet' if not from_port else from_port
        self.units[from_port].add_destination(to_port)

    def add_source(self, from_port, to_port=''):
        """Registers a source unit port 'from_port' with the inlet port of this
        unit, denoted by 'to_port'.

        Args:
            from_port (UnitPort): the port on the source unit
            to_port (str, optional): the name of the port on this unit to
                associate with the destination unit port. If not specified,
                then a default name of 'inlet' is used.

        Returns:
            None
        """
        to_port = 'inlet' if not to_port else to_port
        self.units[to_port].add_source(from_port)

    def add_unit(self, o):
        """
        Adds the unit 'o' to the flowsheet

        Args:
            o (UnitModel): the unit model object

        Returns:
            UnitModel: named unit model object
        """
        self.units[o.unit_name] = o
        setattr(self, o.unit_name, o)
        self.graph.add_node(o.unit_name)
        return o

    def connect(self, from_unit, to_unit, from_port='', to_port='', name=None, doc=None):
        """Defines a connection from one unit to another

        Arguments 'from_port' and 'to_port' are for units like flash tanks
        where different inlet and/or outlet ports have different meanings.

        Args:
            from_unit (UnitModel): source unit for connection stream
            to_unit (UnitModel): destination unit for connection stream
            from_port (str, optional): port name on source unit
            to_port (str, optional): port name on destination unit
            name (str, optional): name for this connector

        Returns:
            None
        """
        if name is None:
            # If name for conection constraint not specified, generate one.
            # This name will take the form:
            # from_unit_source_port__to__to_unit_destination_port
            name = "{}{}__to__{}{}".format(
                from_unit.unit_name,
                '_' + from_port if from_port is not None else '',
                to_unit.unit_name,
                '_' + to_port if to_port is not None else '')

        # Quick fix: for fixed flowsheets, look for from_unit.from_port and
        # to_unit.to_port. If they are found and are both connector objects,
        # then simply link them immediate via a connector.
        tp_con = getattr(to_unit, to_port, None)
        fp_con = getattr(from_unit, from_port, None)
        if isinstance(fp_con, Connector) and isinstance(tp_con, Connector):
            # connector to connector
            if tp_con.port_type != fp_con.port_type:
                # Should imporve message
                raise Exception("Connecting incompatible port types")
            c = Constraint(expr=tp_con == fp_con)
            c.connection = tp_con.port_type
            self.link_map[name] = {
                'from': from_unit, 'to': to_unit,
                'from_port': from_port, 'to_port': to_port, 'doc': doc}
            setattr(self, name, c)
        else:
            # otherwise, register the link in a dictionary
            # link_map is roughly equivalent to John's self.connections
            # Need to decide on a name later
            self.link_map[name] = {
                'from': from_unit, 'to': to_unit, 'from_port': from_port, 'to_port': to_port, 'doc': doc}
            # Register the new connection with both the source and destination
            # ports
            from_unit.add_destination(to_unit.get_in_port(to_port), from_port)
            to_unit.add_source(from_unit.get_out_port(from_port), to_port)

        # Add edge in process map
        self.graph.add_edge(from_unit, to_unit)

    def expand_connectors(self):
        """
        Expands the connectors out to equality constraints. Assuming
        the flowsheet has connectors, this must be done before solving
        the flowsheet model.

        Returns:
            None
        """
        # TODO (from Qi): maybe move build_links() invocation in here?
        xfrm = ConnectorExpander()
        xfrm.apply(instance=self)

    def build_units(self):
        """Constructs the units currently being associated with the flowsheet

        Returns:
            None
        """
        for unit in itervalues(self.units):
            unit.build()
            unit.build_units()
            # Do not expand the links here, because otherwise bad things will
            # happen :(

    def build_links(self):
        """Constructs the links between units as Connector objects

        Returns:
            None
        """
        for name, link in iteritems(self.link_map):
            # Iterate through the dictionary of links
            try:
                # Assume that the link is from a major unit model to another
                # major unit model. Attempt to build the link.
                from_port_name = link['from_port']
                fp_con = link['from'].get_out_port(from_port_name).get_external_link(link)
                to_port_name = link['to_port']
                tp_con = link['to'].get_in_port(to_port_name).get_external_link(link)
                # Need now to check and make sure that the two connectors don't
                # have a mismatch. If they do, then we need to resolve it. The
                # way to resolve this is to apply the following strategy: 1) if
                # the fc vs. F and y are mismatched, and one of them has all of
                # fc, F, and y, then just use the one that's common to both of
                # them. Otherwise, if fc vs. F and y are exclusively
                # mismatched, connect them all and tell the respective ports to
                # turn their linking constraints on. 2) if one port is missing
                # T, P, or vf, then simply remove that variable from the other
                # connector.
                common_elements = set(fp_con.vars) & set(tp_con.vars)
                flows = set(['F', 'y', 'fc'])
                common_flows = common_elements & flows
                if not common_flows:
                    # no common flow elements: choose a scheme, activate linking
                    raise NotImplementedError(
                        'Ability to connect different types of flow elements'
                        ' not supported yet')
                elif common_flows == flows:
                    # all elements are included. May want to select one or the
                    # other to pass through, but for now, just use both.
                    pass
                else:
                    # there are common flows. Remove any extraneous flows.
                    uncommon = set(fp_con.vars) & flows - common_flows
                    for element in uncommon:
                        # remove extraneous flow connector elements
                        del fp_con.vars[element]
                    uncommon = set(tp_con.vars) & flows - common_flows
                    for element in uncommon:
                        # remove extraneous flow connector elements
                        del tp_con.vars[element]
                # if there are any T, P, or vf elements in one connector, but
                # not the other, remove them.
                intrinsic = set(['T', 'P', 'vf'])
                for element in set(fp_con.vars) & intrinsic - common_elements:
                    del fp_con.vars[element]
                for element in set(tp_con.vars) & intrinsic - common_elements:
                    del tp_con.vars[element]
                setattr(
                    self.links, name, pe.Constraint(expr=fp_con == tp_con)
                )
            except KeyError:
                # If the link didn't contain expected 'from' and/or 'to' keys,
                # then it is a different type of link. Assume that the link is
                # from a variable port (Pyomo Connector) to the UnitPort
                # associated with it. Attempt to build the link.
                setattr(self.links, name, pe.Constraint(
                    expr=link['var_port'] == link['unit_port']))
        for unit in itervalues(self.units):
            unit.build_links()

    def apply_OA_strategy(self, oa_ports=False):
        nl_cons = set(self.get_nonlinear_constraints(active=True))
        if len(nl_cons) > 0:
            # only set up OA infrastructure if there are nonlinear constraints
            if oa_ports is False and 'unit_port' in self.unit_type:
                # if we are not applying OA cuts to unit ports and this unit is
                # a unit port, then obviously, we should not apply OA cuts to
                # it
                return
            oa_block = getattr(self, 'oa_cuts', None)
            if oa_block is None:
                # Set up special block for holding OA equations
                oa_block = self.oa_cuts = pe.Block()
            jacs = getattr(self, 'jacs', None)
            if jacs is None:
                jacs = self.jacs = {}
            apply_OA(nl_cons, oa_block, jacs)
            for unit in itervalues(self.units):
                unit.apply_OA_strategy(oa_ports)

    def add_oa_cut(self, iter_num):
        oa_block = getattr(self, 'oa_cuts', None)
        if oa_block is None:
            return  # do not add cuts if OA strategy isn't applied
        exists = getattr(self, 'equip_exists', 1.0)
        if abs(pe.value(exists) - 1.0) <= 1E-6:
            # Unit is active
            add_oa_constraints(
                iter_num, oa_block, self.jacs,
                getattr(self, 'equip_exists', 1.0))
            for unit in itervalues(self.units):
                unit.add_oa_cut(iter_num)

    def activate_oa_cuts(self):
        try:
            self.oa_cuts.activate()
        except AttributeError:
            pass

    def deactivate_oa_cuts(self):
        try:
            self.oa_cuts.deactivate()
        except AttributeError:
            pass

    def build(self):
        raise Exception("ProcessBase build() should be overridden")

    def _reg_initialize(self, name, func, doc=""):
        """
        Register an initialization method for a Pyomo block.  Pyomo
        block initialization methods are used to create a good initial
        starting intial values for variables.
        """
        self._init_methods[name] = _InitMethod(func, doc)

    def display_init_methods(self):
        """
        Print available initialization functions for the Pyomo block and
        their documentation.
        """
        for m in _init_methods:
            print("{0}: {1}".format(m, self._init_methods[m].doc))

    def initialize(self, method, *args, **kwargs):
        """
        Run the specified initialization method on the block.

        The initialization method may work in-place or return a
        dictionary of variable values.  If the initialization method
        returns d, this loads thos values into the Pyomo block.

        Args:
        -------------
        method = mthod key from dict of Pyomo block initialization
                 methods
        *args and **kwargs, any aditional arguments will get passed on
                 to the intialization function.
        """
        self._init_methods[method].func(*args, **kwargs)

    def get_out_port(self, port_name=''):
        """Returns an outlet port object associated with this
        UnitModel

        Args:
            port_name (str, optional): name of the outlet port

        Returns:
            OutPort: outlet port
        """
        if not port_name:
            # if port_name is None, False, or empty, set default value
            port_name = 'outlet'
        return getattr(self.units, port_name)

    def get_in_port(self, port_name=''):
        """Returns an inlet port object associated with this UnitModel

        Args:
            port_name (str, optional): name of the inlet port

        Returns:
            InPort: inlet port
        """
        if not port_name:
            # if port_name is None, False, or empty, set default value
            port_name = 'inlet'
        return getattr(self.units, port_name)

    def solve(self, *args, **kwargs):
        """
        Solve the model.

        Args:
            *args (tuple): Arguments to pass on to the Pyomo solver object
            **kwargs (dict): Arguments to pass on to the Pyomo solver object

        Raises:
            SolverError: if an available solver could not be found
        """
        # Pop the strip bounds argument that is not for Pyomo's solve
        sb = kwargs.pop("strip_bounds", False)
        solver_name = kwargs.pop("using", "local_NLP")
        m = kwargs.pop("model", self)
        try:
            # if the specified solver is in the set of known solvers or solver
            # classes
            solver_spec = self.solvers[solver_name]
            try:
                # assume that we are working with a solver class. Try to find
                # an available solver.
                for sol in solver_spec._solver_list:
                    self.pyomo_solver = SolverFactory(sol)
                    if not self.pyomo_solver.available():
                        self.pyomo_solver = None
                    else:
                        solver_name = sol
                        break
                if self.pyomo_solver is None:
                    # unable to find available solver
                    raise SolverError(
                        'Unable to find available solver among list: {}'
                        .format(solver_spec._solver_list))
                # populate with options
                for key, opt in solver_spec.solver_options(s_name=solver_name):
                    self.pyomo_solver.options[key] = opt
            except AttributeError:
                # It was not a solver class. Just populate with the specific
                # solver options.
                self.pyomo_solver = SolverFactory(solver_name)
                for key, opt in iteritems(solver_spec):
                    self.pyomo_solver.options[key] = opt
        except KeyError:
            self.pyomo_solver = SolverFactory(solver_name)
        if not self.pyomo_solver.available():
            raise SolverError('Solver {} is not available on this system.'
                              .format(solver_name))
        # Set the start time
        self.solve_time = [datetime.now()]
        # Strip the bounds on variables before solving if you want to
        if sb:
            m.strip_bounds()
        self.pyomo_results = self.pyomo_solver.solve(m, *args, **kwargs)
        # Put bounds back on if they were stripped off.
        if sb:
            m.restore_bounds()
        # Record solver end time
        self.solve_time.append(datetime.now())
        # Return the results of the solve
        # TODO: remove temporary fix for Pyomo bug <JCE>
#        try:
#            del os.environ['AMPLFUNC']
#        except KeyError:
#            pass
        return self.pyomo_results

    def display_results(self):
        """
        Show solve result info
        """
        print(self.pyomo_results)

    def strip_bounds(self):
        """
        Print variables that are near bounds.  Marks variables close
        to bounds.
        """
        self.bstore = j2.save_json(self, wts=j2.StoreSpecBounds())
        blk = self
        for o in blk.component_objects(descend_into=True):
            if isinstance(o, Var):
                for i, v in o.iteritems():
                    v.setlb(None)
                    v.setub(None)

    def restore_bounds(self):
        j2.load_json(self, sd=self.bstore, wts=j2.StoreSpecBounds())

    def large_residuals(self, tol=1e-5):
        for o in self.component_objects(descend_into=True):
            if isinstance(o, pyomo.core.base.constraint.IndexedConstraint):
                for i, c in o.iteritems():
                    if c.lower - c.body() > 1e-5 and c.active:
                        print("{} {}".format(c.name, c.lower - c.body()))
            elif isinstance(o, Constraint):
                c = o
                try:
                    if c.lower - c.body() > 1e-5 and c.active:
                        print("{} {}".format(c.name, c.lower - c.body()))
                except:
                    pass

    def unit_ports(self):
        """Returns a generator over the unit ports in the unit model. Does not
        attempt to recurse through sub-units for their unit ports.

        Returns:
            generator: generator over unit ports in the unit model
        """
        return (unit for unit in itervalues(self.units)
                if 'unit_port' in unit.unit_type)

    def all_units(self):
        """Generator function yielding all units and their respective sub-units

        Returns:
            UnitModel: unit model
        """
        for unit in itervalues(self.units):
            yield unit
            for subunit in unit.all_units():
                yield subunit

    def fetch_limits(self):
        """Sets the min and max values for certain properties for a unit based
        upon the values of its containing (parent) unit

        Returns:
            None
        """
        if self.parent_unit is None:
            return
        self.max_flow = self.parent_unit.max_flow
        self.min_T = self.parent_unit.min_T
        self.max_T = self.parent_unit.max_T
        self.max_P = self.parent_unit.max_P

    @property
    def max_flow(self):
        """
        Sets a default upper bound on the value of all flow variables in this
        unit as well as the units it contains. No guarantee is made that its
        descendants will adhere to this maximum, but by default, this will be
        true.

        This should be dictated, when possible, by physical insight; for
        example, by practical limits on pipe size. Tighter variable bounds can
        be set on a per-unit basis by other mechanisms.

        Note: this limit is also placed on 'total_flow' variables.

        Returns:
            float: maximum flowrate [units TODO]
        """
        try:
            return self._max_flow
        except AttributeError:
            return None

    @max_flow.setter
    def max_flow(self, value):
        self._max_flow = pe.value(value)
        for unit in itervalues(self.units):
            unit.max_flow = value

    @property
    def min_T(self):
        """Minimum temperature for this unit. No guarantee is made that its
        descendants will adhere to this maximum, but by default, this will be
        true.

        Returns:
            float: minimum temperature [units TODO]
        """
        try:
            return self._min_T
        except AttributeError:
            return None

    @min_T.setter
    def min_T(self, value):
        self._min_T = pe.value(value)
        for unit in itervalues(self.units):
            unit.min_T = value

    @property
    def max_T(self):
        """Maximum temperature for this unit. No guarantee is made that its
        descendants will adhere to this maximum, but by default, this will be
        true.

        Returns:
            float: maximum temperature [units TODO]
        """
        try:
            return self._max_T
        except AttributeError:
            return None

    @max_T.setter
    def max_T(self, value):
        self._max_T = pe.value(value)
        for unit in itervalues(self.units):
            unit.max_T = value

    @property
    def max_P(self):
        """Maximum pressure for this unit. No guarantee is made that its
        descendants will adhere to this maximum, but by default, this will be
        true.

        Returns:
            float: maximum pressure [units TODO]
        """
        try:
            return self._max_P
        except AttributeError:
            return None

    @max_P.setter
    def max_P(self, value):
        self._max_P = pe.value(value)
        for unit in itervalues(self.units):
            unit.max_P = value
