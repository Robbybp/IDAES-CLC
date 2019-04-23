"""
Base clase for unit models
"""
from __future__ import division
from __future__ import print_function

from .process_base import ProcessBase
import pyomo.environ as pe
from six import iteritems, itervalues
from .util.var import tighten_var_bound
from .util.cut_gen import clone_block, get_sum_sq_diff, count_vars
from .util.var import lb, ub

__author__ = "John Eslick <john.eslick@netl.doe.gov>, "\
             "Qi Chen <qichen@andrew.cmu.edu>"
__version__ = "1.0.0"

__all__ = ['UnitModel']


class UnitModel(ProcessBase):
    """
    This is the class for process unit operations models. These are models that
    would generally appear in a process flowsheet or superstructure. Aside from
    a default model initialization function, the only thing here that is not in
    the ProcessBase class is the methods to make standard connectors.

    Attributes:
        block_bounds (dict): Dictionary of variable bounds that should be satisfied when the unit is active. See __init__.
        comps (Set): Pyomo Set of components in the process
        equip (Block): Pyomo Block containing ordinary unit model constraints
        jacs (dict): Dictionary of jacobians
        lin_cuts (Block): Pyomo Block containing generated rigorous linear constraints
        max_flow (float): maximum value for any flow variables
        max_P (float): maximum value for any pressure variables
        max_T (float): maximum value for any temperature variables
        min_T (float): minimum value for any temperature variables
        parent_unit (Wrapper): Wrapper around the UnitModel that this unit is
            contained within
        unit_name (str): Name of the process unit
        unit_type (set): set of tags describing the unit
        var_bounds (dict): Dictionary of variable bounds that should always be
            active. See __init__.
    """

    def __init__(self, *args, **kwargs):
        """Initializes the UnitModel, which represents a flowsheet unit

        Kwargs:
            type (list): List of strings describing the unit type
            bounds (dict): dictionary mapping _VarData local_name to
            (lower bound, upper bound) tuples for variable bounds that should
            always be satisfied
            block_bounds (dict): dictionary mapping _VarData local_name to
            (lower bound, upper bound) tuples for variable bounds that should
            be satisfied when the unit is active
        """
        kwargs.setdefault('name', 'Unnamed_Unit_Model')

        # variables that are temporarily fixed during LOA
        self._tmp_fixed = set()
        # nonlinear constraints that are temporarily deactivated during LOA
        self._tmp_nonlinear_deactivated = set()
        # flow variables upon which a temporary lower bound is enforced
        self._tmp_min_flow = set()

        self.block_bounds = kwargs.pop('block_bounds', {})
        # convert tuples to dictionaries
        for k, v in iteritems(self.block_bounds):
            self.block_bounds[k] = {'lb': v[0], 'ub': v[1]}

        self.var_bounds = kwargs.pop('bounds', {})
        self.jacs = {}
        self._fathomed = False

        super(UnitModel, self).__init__(*args, **kwargs)

    def build(self):
        """Constructs the UnitModel

        Performs the following in sequential order:
         1) constructs an 'equip' sub-block for unit constraints
         2) constructs an 'lin_cuts' sub-block for rigorous linear constraints
         3) calls the '_build_unit_sets()' function
         4) calls the '_build_unit_params()' function
         5) calls the '_build_unit_vars()' function
         6) calls the '_build_unit_constraints()' function
         7) calls the '_tighten_bounds()' function
         8) calls the '_enforce_block_bounds()' function
         9) calls the '_build_unit_ports()' function

        The '_build_unit_[element]' functions above may be overriden by
        subclasses in order to perform additional construction tasks.

        Returns:
            None
        """
        if self._built:
            return
        self.equip = pe.Block()
        self.equip.block_bounds = pe.ConstraintList()
        # self.oa_cuts = pe.Block()  # constructed if necessary in
        # apply_OA_strategy() function
        self.lin_cuts = pe.Block()
        self.lin_cuts.self_proj_cuts = pe.ConstraintList()
        self._build_unit_sets()
        self._build_unit_params()
        self._build_unit_vars()
        self._build_unit_constraints()
        self._tighten_bounds()
        # TODO validate that block_bounds passed in refer to actual variables
        self._enforce_block_bounds()
        # give warning otherwise
        self._build_unit_ports()
        self._build_unit_initialization()
        self._built = True

    @property
    def unit_name(self):
        return self._unit_name

    @unit_name.setter
    def unit_name(self, new_name):
        self._unit_name = new_name

    def get_nonlinear_constraints(self, **kwargs):
        """Generator function over the nonlinear constraints in the unit model.

        Defaults for component_data_objects are used, unless passed in as
        keyword arguments otherwise. As of this writing, this function will
        yield inactive constraints but will not recurse into sub-blocks.

        Returns:
            _ConstraintData: generator over _ConstraintData objects

        Raises:
            ValueError: if polynomial degree is not positive
        """
        for constr in self.equip.component_data_objects(
                ctype=pe.Constraint, **kwargs):
            deg = constr.body.polynomial_degree()
            if deg is None or deg > 1:
                yield constr
            elif deg == 0 or deg == 1:
                # constraint is linear
                pass
            else:
                raise ValueError("Unexpected polynomial degree: " + deg)

    def reactivate_nonlinear_constraints(self):
        for constrID in self._tmp_nonlinear_deactivated:
            constrID.find_component_on(self.model()).activate()
        self._tmp_nonlinear_deactivated.clear()

    def deactivate_nonlinear_constraints(self):
        for constr in self.get_nonlinear_constraints(active=True):
            constr.deactivate()
            self._tmp_nonlinear_deactivated.add(pe.ComponentUID(constr))

    def _build_unit_sets(self):
        """Subclasses should override this function to construct necessary Set
        objects
        """
        pass

    def _build_unit_params(self):
        """Subclasses should override this function to construct necessary
        Param objects
        """
        pass

    def _build_unit_vars(self):
        """Subclasses should override this function to construct necessary Var
        objects
        """
        pass

    def _build_unit_constraints(self):
        """Subclasses should override this function to construct necessary
        Constraint objects
        """
        pass

    def _build_unit_ports(self):
        """Subclasses should override this function to construct necessary Port
        objects
        """
        pass

    def _build_unit_initialization(self):
        """Subclasses may override this function in order to perform variable
        initialization tasks immediately post-build.

        Returns:
            None
        """
        pass

    def _tighten_bounds(self):
        """Tightens the variable bounds based upon information passed to
        the UnitModel at initialization in the 'bounds' dictionary

        Returns:
            None
        """
        if not self.var_bounds:
            return  # if bounds dictionary is empty, skip all this.
        bound_vars = set()
        # Set prespecified bounds passed in during initialization
        for vardata in self.component_data_objects(ctype=pe.Var):
            var_item_name = vardata.local_name
            var_name = vardata.parent_component().local_name
            if var_item_name in self.var_bounds:
                bound_vars.add(var_item_name)
                bound_vars.add(var_name)
                tighten_var_bound(vardata, self.var_bounds[var_item_name])
            elif var_name in self.var_bounds:
                bound_vars.add(var_name)
                tighten_var_bound(vardata, self.var_bounds[var_name])
            else:
                continue
        missed_bounds = set(self.var_bounds) - bound_vars
        for bnd in missed_bounds:
            print('Warning: variable {} was not found in {}'
                  .format(bnd, self.local_name))

    def _enforce_block_bounds(self):
        """Generates a ConstraintList enforcing the active disjunction variable
        bounds explicitly

        Deletes the ConstraintList 'self.equip.block_bounds' if it
        already exists and populates a new ConstraintList with the current
        values in the 'self.block_bounds' dictionary.

        Returns:
            None
        """
        if not self.block_bounds:
            return  # if block_bounds dictionary is empty, skip.
        # if this equipment does not have an equip_exists variable, skip.
        exists = getattr(self, 'equip_exists', None)
        if exists is None:
            print('Warning: block bounds are defined for {}, but the unit '
                  'does not support them.'.format(self.local_name))
            return
        # delete the existing ConstraintList, if it exists.
        try:
            self.equip.block_bounds.clear()
        except:
            self.equip.block_bounds = pe.ConstraintList()
        bb = self.equip.block_bounds
        # Iterate through all vardata objects, adding requisite bounds to the
        # ConstraintList
        for vardata in self.component_data_objects(ctype=pe.Var):
            if vardata.local_name in self.block_bounds:
                # check if block_LB for this vardata is tighter than the LB
                block_LB = self.block_bounds[vardata.local_name]['lb']
                if block_LB is not None:
                    if vardata.lb is None:
                        print('Warning: cannot generate block bound for {}'
                              ' in {} due to missing LB'
                              .format(vardata.local_name, self.local_name))
                    elif block_LB > vardata.lb:
                        bb.add(expr=vardata >= block_LB *
                               exists + vardata.lb * (1 - exists))
                # check if block_UB for this vardata is tighter than the UB
                block_UB = self.block_bounds[vardata.local_name]['ub']
                if block_UB is not None:
                    if vardata.ub is None:
                        print('Warning: cannot generate block bound for {}'
                              ' in {} due to missing UB'
                              .format(vardata.local_name, self.local_name))
                    elif block_UB < vardata.ub:
                        bb.add(expr=vardata <= block_UB *
                               exists + vardata.ub * (1 - exists))

    def apply_linear_relaxations(self, nsegs=1, recurse=True):
        """
        While it would be nice to someday have the logic for applying linear
        relaxations in this method, at present subclasses need to override this
        to provide their own implementation.

        Subclasses also need to take care not to forget to pass on the call to
        their respective children units.
        """
        self._lin_cuts_nsegs = nsegs
        if recurse:
            for unit in itervalues(self.units):
                unit.apply_linear_relaxations(nsegs=nsegs)

    def reconstruct_envelopes(self):
        self.del_component('lin_cuts')
        self.lin_cuts = pe.Block()
        for unit in itervalues(self.units):
            unit.reconstruct_envelopes()
        self.apply_linear_relaxations(nsegs=self._lin_cuts_nsegs, recurse=False)

    def apply_MIP(self):
        """Prepares the UnitModel to be solved as a MILP

        Unfixes the 'equip_exists' variable.
        Deactivates the nonlinear constraints associated with the unit and
        activates the linearization constraints.

        Returns:
            None
        """
        if not self._fathomed:
            self.equip.activate()
            self.unfix_flows()
            try:
                self.equip_exists.unfix()
            except AttributeError:
                pass
            self.lin_cuts.activate()
            self.activate_oa_cuts()
            self.deactivate_nonlinear_constraints()
            for unit in itervalues(self.units):
                unit.apply_MIP()
        else:
            try:
                self.equip_exists.fix(0)
            except AttributeError:
                pass
            self.lin_cuts.deactivate()
            self.deactivate_oa_cuts()
            self.equip.deactivate()
            self.fix_flows()

    def apply_NLP(self):
        """Prepares the UnitModel to be solved as an NLP
        """
        exists = getattr(self, 'equip_exists', 1.0)
        try:
            exists.fix()
        except AttributeError:
            pass
        self.lin_cuts.deactivate()
        self.deactivate_oa_cuts()
        for unit in itervalues(self.units):
            unit.apply_NLP()
        if abs(pe.value(exists) - 1.0) <= 1E-6:
            self.equip.activate()
            self.reactivate_nonlinear_constraints()
            self.unfix_flows()
        else:
            self.equip.deactivate()
            self.fix_flows()

    def fix_flows(self):
        for var in self._get_flow_vars():
            if var.fixed:
                pass
            else:
                self._tmp_fixed.add(pe.ComponentUID(var))
                var.fix(0)

    def unfix_flows(self):
        for varUID in self._tmp_fixed:
            var = varUID.find_component_on(self.model())
            var.unfix()
        self._tmp_fixed.clear()

    def introspect_flows(self):
        """This function examines the flow values of attached connections for
        those that are deactivated, and fixes the corresponding flow-related
        variables to zero.

        This function is really meant mostly for ports, so by default, just
        pass on the call to its child units.

        Returns:
            None
        """
        for unit in itervalues(self.units):
            unit.introspect_flows()

    def deactivate_trivial_constraints(self):
        """Find and deactivates trivial constraints: those that are const =
        const, or const < const. Verifies that they evaluate to True first.

        Stores deactivated constraint UIDs in the set
        '_tmp_trivial_deactivated'

        Returns:
            None
        """
        if not hasattr(self, '_tmp_trivial_deactivated'):
            self._tmp_trivial_deactivated = set()
        for condata in self.component_data_objects(ctype=pe.Constraint, active=True, descend_into=True):
            # if the constraint is trivial, deactivate it.
            if condata.body.is_fixed():
                if not (condata.upper is None or condata.upper.is_fixed()) or not (condata.lower is None or condata.lower.is_fixed()):
                    # that's weird, why would these not be fixed?
                    raise NotImplementedError('Non-fixed upper or lower bound on constraint.')
                if pe.value(condata.body) != pe.value(condata.upper) or pe.value(condata.body) != pe.value(condata.lower):
                    # inconsistent or inequality. Check if inequality.
                    if condata.upper is None:
                        # Make sure that body >= lower
                        if pe.value(condata.body) < pe.value(condata.lower):
                            raise ValueError('Infeasible constraint: ' + condata.name)
                    elif condata.lower is None:
                        # Make sure that body <= upper
                        if pe.value(condata.body) > pe.value(condata.upper):
                            raise ValueError('Infeasible constraint: ' + condata.name)
                    else:
                        raise NotImplementedError('Infeasible constraint: ' + condata.name)
                # otherwise, it's fine
                condata.deactivate()
                self._tmp_trivial_deactivated.add(pe.ComponentUID(condata))

    def reset_trivial_constraints(self):
        """Resets the constraints deactivated as trivial earlier.

        Returns:
            None
        """
        for varUID in self._tmp_trivial_deactivated:
            var = varUID.find_component_on(self.model())
            var.activate()
        self._tmp_trivial_deactivated.clear()

    def reset_introspect_fixed(self):
        """Resets variables fixed during introspection

        Returns:
            None
        """
        for unit in itervalues(self.units):
            unit.reset_introspect_fixed()

    def set_min_flows(self):
        min_flow = 1E-8
        for var in self._get_flow_vars():
            if not var.fixed and lb(var) == 0 and ub(var) >= min_flow:
                self._tmp_min_flow.add(pe.ComponentUID(var))
                var.setlb(min_flow)
        for unit in itervalues(self.units):
            unit.set_min_flows()

    def reset_min_flows(self):
        for varUID in self._tmp_min_flow:
            var = varUID.find_component_on(self.model())
            var.setlb(0)
        self._tmp_min_flow.clear()
        for unit in itervalues(self.units):
            unit.reset_min_flows()

    def get_flow_vars(self):
        raise NotImplementedError(
            'Subclass {} needs to override this generator method to yield '
            'the flow variables.'.format(type(self).__name__))

    def _get_flow_vars(self):
        """get_flow_vars may return indexed vars. This simply turns them all
        into VarData

        Returns:
            _VarData: generator of flow variables
        """
        for var in self.get_flow_vars():
            try:
                # iterate through the variable index
                for indx in var:
                    yield var[indx]
            except TypeError:
                # variable is not iterable. Must be a _VarData
                yield var

    def get_slack_variables(self):
        try:
            for v in self.oa_cuts.component_data_objects(ctype=pe.Var):
                if v.local_name.startswith('_slack_'):
                    yield v
        except AttributeError:
            pass
        for unit in itervalues(self.units):
            for v in unit.get_slack_variables():
                yield v

    def display_flows(self):
        """Displays component flow variables associated with the UnitModel

        Returns:
            None
        """
        for var in self.component_objects(ctype=pe.Var):
            if var.local_name.startswith('flow_') or var.local_name == 'flow' \
                    or var.local_name.startswith('fc_'):
                var.display()

    def display_conc(self):
        """Displays component concentrations associated with the UnitModel

        Returns:
            None
        """
        for var in self.component_objects(ctype=pe.Var):
            if var.local_name.startswith('conc_') or var.local_name == 'conc':
                var.display()

    def display_total_flows(self):
        """Displays total flow variables associated with the UnitModel

        Returns:
            None
        """
        for var in self.component_objects(ctype=pe.Var):
            if var.local_name.startswith('total_flow_') or \
                    var.local_name == 'total_flow':
                var.display()

    def display_T(self):
        """Displays temperature variables associated with the UnitModel

        Returns:
            None
        """
        for var in self.component_objects(ctype=pe.Var):
            if var.local_name.startswith('T_') or \
                    var.local_name in ('T', 'Tin', 'Tout'):
                var.display()

    def display_P(self):
        """Displays pressure variables associated with the UnitModel

        Returns:
            None
        """
        for var in self.component_objects(ctype=pe.Var):
            if var.local_name.startswith('P_') or \
                    var.local_name in ('P', 'Pin', 'Pout'):
                var.display()

    def display_variables(self, simple=False, descend_into=True):
        """Displays all variables associated with the UnitModel

        Args:
            simple (bool, optional): Print a simplified version showing only
                variable values.

        Returns:
            None
        """
        if not simple:
            for var in self.component_objects(ctype=pe.Var, descend_into=descend_into):
                var.display()
        else:
            for vardata in self.component_data_objects(ctype=pe.Var, descend_into=False):
                print("{}: {}".format(vardata.local_name, vardata.value))

    def display_base_variables(self, simple=False):
        """Displays all variables at the root block of the UnitModel

        Args:
            simple (bool, optional): Print a simplified version showing only
                variable values.

        Returns:
            None
        """
        self.display_variables(simple=simple, descend_into=False)

    def generate_cut_gen_problem(self):
        if not hasattr(self, 'equip_exists'):
            pass
        elif abs(self.equip_exists.value) <= 1E-3:
            # Do not solve cut generation problem for inactive units
            return None

        if not hasattr(self, 'apply_NLP'):
            return None

        self.apply_NLP()
        num_vars = count_vars(self)
        if num_vars < 1:
            return None

        b = pe.ConcreteModel(name=self.local_name)
        # self.cg_prob = b
        var_set = b.var_set = pe.RangeSet(num_vars)
        b.lbda = pe.Var(b.var_set, domain=pe.NonNegativeReals,
                        bounds=(0, 1), initialize=0)
        b.lbda_sum = pe.Constraint(
            expr=sum(b.lbda[vs] for vs in b.var_set) == 1)

        clone_block(self, b, var_set, b.lbda)

        sum_sq_diff = get_sum_sq_diff(self, b)
        b.obj = pe.Objective(expr=sum_sq_diff, sense=pe.minimize)
        return b

    def apply_self_proj_cut(self, cg_prob):
        # Note: this is still buggy
        if pe.value(cg_prob.obj.expr) >= 1E-3:
            # generate cut
            print(self.name, pe.value(cg_prob.obj.expr))
            # self.cg_prob = cg_prob
            print('Adding self-projection cut for ' + self.local_name)
            self.lin_cuts.self_proj_cuts.add(
                2 * sum(
                    (getattr(cg_prob, '_{}_clone'.format(
                        var.parent_component().local_name))[var.index()].value -
                     var.value
                     ) *
                    (var -
                     getattr(cg_prob, '_{}_clone'.format(
                         var.parent_component().local_name))[var.index()].value
                     )
                    for var in self.component_data_objects(
                        ctype=pe.Var, active=True) if not var.fixed
                ) >= 0
            )
        else:
            print('Self-projection cut for ' +
                  self.local_name + ' redundant.')

    def fathom(self):
        """Indicates that this unit cannot be active in the optimal
        superstructure configuration. Fix the unit to inactive and do not allow
        it to be active in the MIP master problem.

        Returns:
            None
        """
        self._fathomed = True
        try:
            self.equip_exists.fix(0)
        except AttributeError:
            pass
        for unit in itervalues(self.units):
            unit.fathom()

    def fetch_limits(self):
        """Sets the min and max values for certain properties for a unit based
        upon the values of its containing (parent) unit

        Returns:
            None
        """
        self.max_flow = self.parent_unit.max_flow
        self.min_T = self.parent_unit.min_T
        self.max_T = self.parent_unit.max_T
        self.max_P = self.parent_unit.max_P
