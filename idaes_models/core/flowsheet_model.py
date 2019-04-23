"""
This is a class for flowsheet (AKA process or main) models.
"""
from __future__ import division
from __future__ import print_function

import pyomo.environ as pe
from pyomo.opt import TerminationCondition, SolutionStatus, SolverStatus
# from pyomo.opt.parallel import SolverManagerFactory
from six import itervalues, iterkeys
from idaes_models.core.util.misc import doNothing, get_time, rround
from bunch import Bunch

from idaes_models.core.process_base import ProcessBase
from idaes_models.core.plugins.propagate_fixed import propagate_var_fix, reset_propagated_var_fix

import random

# Some more inforation about this module
__author__ = "John Eslick <john.eslick@netl.doe.gov>, "\
             "Qi Chen <qichen@andrew.cmu.edu>"
__version__ = "1.0.0"

__all__ = ['FlowsheetModel']


class FlowsheetModel(ProcessBase):
    """This object stores and manages the process synthesis flowsheet model.

    Attributes:
        comps (list): component list
        links (dict): dictionary of connections between process units
        logic (dict): dictionary of logical relations between process units
        name (string): name of the flowsheet
        oa_iter (int): counter for number of Outer Approximation iterations
        results (SolverResults): holds the most recent solver results object
        solver_time (dict): holds a dictionary of floats for the number of
            seconds that different subproblems take
        units (dict): dictionary holding all process units
        use_UB_constr (bool): flag for whether to enforce the global upper
            bound as a constraint
    """

    def __init__(self, *args, **kwargs):
        """Initializes a new FlowsheetModel

        Kwargs:
            name (str): Name of the flowsheet
        """
        kwargs.setdefault('name', 'Unnamed_Flowsheet')
        kwargs.setdefault('parent', None)
        super(FlowsheetModel, self).__init__(*args, **kwargs)

    def __init2__(self, *args, **kwargs):
        super(FlowsheetModel, self).__init2__()
        self.int_cuts = pe.ConstraintList()
        self.tmp_int_cuts = pe.ConstraintList()
        self._global_LB =\
            pe.Param(initialize=float('-inf'), mutable=True)
        self._local_LB = float('-inf')
        self._global_UB =\
            pe.Param(initialize=float('inf'), mutable=True)
        self.use_UB_constr = False
        self.solver_time = Bunch(mip=0, local_nlp=0, global_nlp=0,
                                 self_proj_cut=0, bounding=0, global_minlp=0)
        self.dual = pe.Suffix(direction=pe.Suffix.IMPORT_EXPORT)
        self.dual.deactivate()
        self.oa_iter = 0
        self._best_found_result = None
        self.logic = Bunch()
        self.logic.setdefault('one_of', set())
        self.logic.setdefault('at_least_one_of', set())
        self.logic.setdefault('all_of', set())

    def build(self):
        pass

    @property
    def global_LB(self):
        """The global lower bound on the objective value.

        When maximization is desired, then the negative of the original
        objective value will be considered.

        If no lower bound is available, then a value of float('-inf') will be
        returned.

        Returns:
            float: global lower bound value
        """
        return pe.value(self._global_LB)

    @property
    def local_LB(self):
        """The non-rigorous lower bound on the objective value.

        Returns:
            float: non-rigorous lower bound value
        """
        return self._local_LB

    @property
    def global_UB(self):
        """Gives the global upper bound on the objective value.

        This value corresponds to the objective value of the best found
        feasible flowsheet configuration. If no such configuration exists,
        then a value of float('inf') will be returned.

        Returns:
            float: global upper bound value
        """
        return pe.value(self._global_UB)

    @property
    def total_solver_time(self):
        """The total solution time among tracked solvers

        Returns:
            float: cumulative CPU time reported by various solvers
        """
        return sum(time for time in itervalues(self.solver_time))

    @property
    def comp_set(self):
        """The chemical components present in the flowsheet.

        Returns:
            None
        """
        return self._comps

    @comp_set.setter
    def comp_set(self, comps):
        """Sets the chemical components present in the flowsheet

        Args:
            comps (list): a list of strings representing chemical components

        Returns:
            None
        """
        self._comps = set(comps)

    def _activate_standard_objective(self):
        """Activates the standard objective.

        Subclasses should override this function to activate the appropriate
        objective function and deactivate all other objective functions.

        Raises:
            NotImplementedError: if not overriden by subclass
        """
        raise NotImplementedError()

    def _activate_penalized_oa_objective(self):
        """Activates the penalized outer approximation objective.

        Subclasses should override this function to activate the appropriate
        objective function and deactivate all other objective functions.

        Raises:
            NotImplementedError: if not overriden by subclass
        """
        raise NotImplementedError()

    def _apply_UB_constr(self):
        """Applies an upper bound constraint if appropriate

        Checks if an upper bound constraint should be imposed upon the
        objective and if so, does that. This is primarily intended as a way to
        "shortcut" global NLP subproblems where the subproblem LB is greater
        than the global UB.

        The current implementation deletes the old constraint, if it exists,
        and replaces it with an updated one.
        """
        if self.use_UB_constr and self.global_UB < float('inf'):
            self.del_component('global_UB_constr')
            self.global_UB_constr = pe.Constraint(
                expr=self.obj.expr <= self.global_UB)

    def solve_global_MINLP(self):
        """Solves the global full-space MINLP problem directly

        First activates the standard objective, then calls the global_MINLP
        solver to solve the problem. By default, keepfiles and tee are both
        set to True.

        Returns:
            None
        """
        # Cost calculation
        self._activate_standard_objective()
        self.pyomo_results = self.solve(
            using='global_MINLP', tee=True, keepfiles=True)
        self.solver_time.global_minlp += self.pyomo_results.solver.time
        self.solutions.store_to(self.pyomo_results)
        # Assumes that this function is not run multiple times
        self._best_found_result = self.pyomo_results
        self._global_UB = min([self.global_UB, pe.value(
            self.pyomo_results.problem.upper_bound)])
        self._global_LB = max([self.global_LB, pe.value(
            self.pyomo_results.problem.lower_bound)])

    def solve_local_NLP(self, tee=False, keepfiles=False, add_min_flows=True, iterations=5):
        """Solves the nonlinear program for a fixed flowsheet configuration

        This NLP is used as part of the LOA algorithm to assess a given
        flowsheet configuration.

        The solver in this case is a local optimization solver, so it is not
        guaranteed to find the best possible solution for the configuration.
        However, any feasible solution found represents a potential upper bound
        on the objective value, so it is updated correspondingly.

        If a feasible solution is found, its values are loaded back into the
        model objects.

        Args:
            tee (bool, optional): flag for verbose solver output
            keepfiles (bool, optional): flag to keep temporary solver files

        Returns:
            bool: True if feasible solution found. False otherwise.
        """
        # Cost calculation
        self._activate_standard_objective()
        for o in itervalues(self.units):
            o.apply_NLP()
        self.dual.activate()
        propagate_var_fix(self, tmp=True)
        for o in itervalues(self.units):
            o.introspect_flows()
        for o in itervalues(self.units):
            o.deactivate_trivial_constraints()
        for o in itervalues(self.units):
            if add_min_flows:
                o.set_min_flows()
        # TODO remove skip_trivial_constraints when Conopt empty rows
        # interface fixed
        results = self.pyomo_results = self.solve(
            using='local_NLP', tee=tee, load_solutions=False, skip_trivial_constraints=True, keepfiles=keepfiles)
        reset_propagated_var_fix(self)
        for o in itervalues(self.units):
            o.reset_introspect_fixed()
            o.reset_trivial_constraints()
        if add_min_flows:
            for o in itervalues(self.units):
                o.reset_min_flows()
        self.solver_time.local_nlp += results.solver.time
        print('Local NLP took ' + str(round(results.solver.time, 2)) + 's')
        if results.solver.status is SolverStatus.ok and \
                results.solver.termination_condition is \
                TerminationCondition.optimal:
            self.solutions.load_from(results)
            if pe.value(self.obj.expr) < self.global_UB:
                self.solutions.store_to(self.pyomo_results)
                self._best_found_result = self.pyomo_results
            self._global_UB = min([self.global_UB, pe.value(self.obj.expr)])
            print('NLP subproblem found local solution of ' +
                  str(pe.value(self.obj.expr)))
            self.dual.deactivate()
            return True
        else:
            # To recover dual information from an infeasible run, uncomment the
            # following two lines.
            # self.solutions.load_from(results)
            # self.dual.display()
            # import sys
            # sys.exit()
            self.dual.deactivate()
            print('NLP subproblem did not converge to local solution.')
            completed = False
            if(iterations > 1):
                self.reinitialize_local_NLP("rand_guess_and_bound")
                completed = self.solve_local_NLP(
                    tee=tee, keepfiles=keepfiles, add_min_flows=add_min_flows, iterations=iterations - 1)
            return completed

    def reinitialize_local_NLP(self,strategy = "rand_distributed"):
        """ reinitializes local NLP in case of non convergence

            If local NLP does not find a solution with the given
            initialization, use the given strategy to find and set
            a new starting value for all variables

            Args:
                strategy(string, optional) : determines the
                    strateguy used for reinitialization.

            Strategies:
                midpoint returns the arithmetic mean of the
                upper and lower bound

                random returns a random point between the
                upper and lower bound (note: can equal lower
                bound

                midpoint_guess_and_bound returns the midpoint
                between the initial guess and the farthest bound

                rand_guess_and_bound returns a random point
                between the initial guess and the farthest
                bound
        """

        def midpoint(val, lb, ub):
            return (lb + ub) // 2

        def rand(val, lb, ub):
            return (ub - lb) * random.random() + lb

        def midpoint_guess_and_bound(val, lb, ub):
            bound = ub if ((ub - val) >= (val - lb)) else lb
            return (bound + val) // 2

        def rand_guess_and_bound(val, lb, ub):
            bound = ub if ((ub - val) >= (val - lb)) else lb
            return (abs(bound - val) * random.random()) + min(bound, val)

        def rand_distributed(val, lb, ub, divisions = 9):
            linspace = []
            spacing = (ub-lb)/divisions
            for i in range(divisions-1):
                linspace[i] = spacing * (i+1)
            return random.choice(linspace)


        #iterate thorugh all units in flowsheet
        print("attempting reinitialization on all variables")
        for o in itervalues(self.units):
            #o.display()
            #o.solve()
            #print("AFTER SOLVING")
            #o.display();
            #iterate through all variables in that uni
            for var in self.component_data_objects(ctype=pe.Var, descend_into=True):
                if not var.is_fixed() and not var.is_binary and not var.is_integer \
                        and not (val == None or lb == None or  ub == None or strategy == None):
                    val = value(var)
                    lb  = var.lb
                    ub  = var.ub
                    #apply strategy to bounds/variable
                    strategies = {"midpoint"    : midpoint,
                                  "random"      : random,
                                  "midpoint_guess_and_bound" : midpoint_guess_and_bound,
                                  "rand_guess_and_bound" : rand_guess_and_bound,
                                  "rand_distributed" : rand_distributed
                                  }

                    new = strategies[strategy](val, lb, ub)
                    var = new

    def solve_MIP(self, tee=False):
        """Solves the full space mixed-integer linear approximation

        The MIP is used as part of the GLOA algorithm to determine a new
        configuration for further assessment. This function assumes that the
        linearizations applied to form the linear approximation can be
        regarded as rigorous and therefore that the objective value from
        solution of the MIP is a global lower bound, which is updated based on
        the solution.

        If the MIP is infeasible, then the lower bound is set to a value of
        float('inf').

        Args:
            tee (bool, optional): flag for verbose solver output

        Returns:
            bool: True if feasible solution found. False otherwise.
        """
        self._activate_standard_objective()
        for o in itervalues(self.units):
            getattr(o, 'apply_MIP', doNothing)()
        results = self.pyomo_results = \
            self.solve(using='global_MIP', tee=tee)
        self.solver_time.mip += get_time(results)
        print('MIP took ' + str(round(get_time(results), 2)) + 's')
        if results.solver.termination_condition is\
                TerminationCondition.optimal:
            self._global_LB = max([self.global_LB,
                                   pe.value(self.obj.expr)])
            print('MIP Master: ' + (str(round(pe.value(self.obj.expr))) if results.solver.termination_condition is TerminationCondition.optimal else 'infeasible') +
                  ' L: ' + str(round(self.global_LB)) + ' U: ' + str(round(self.global_UB)))
            print('MIP proposes a configuration of:')
            self.print_active_units()
            return True
        else:
            print('MIP infeasible')
            self._global_LB = float('inf')
            return False

    def solve_local_MIP(self, tee=False):
        """Solves the non-rigorous lower-bounding MIP problem

        The MIP is used as part of the LOA algorithm to determine a new
        configuration for further assessment. The lower bound provided by
        this function is not rigorous because of OA cuts on potentially
        nonconvex functions. For the OA cuts, the equality relaxation with
        an augmented penalty function method is applied.

        Args:
            tee (bool, optional): flag for verbose solver output

        Returns:
            bool: True if feasible solution found. False otherwise.
        """
        self._activate_standard_objective()
        self._activate_penalized_oa_objective()
        for o in itervalues(self.units):
            getattr(o, 'apply_MIP', doNothing)()
        results = self.pyomo_results = \
            self.solve(using='local_MIP', tee=tee)
        self.solver_time.mip += get_time(results)
        print('MIP took ' + str(round(get_time(results), 2)) + 's')
        if results.solver.termination_condition is\
                TerminationCondition.optimal:
            self._local_LB = max([self.local_LB,
                                  self.global_LB,
                                  pe.value(self.oa_obj.expr)])
            print('MIP Master: ' + (str(round(pe.value(self.oa_obj.expr))) if results.solver.termination_condition is TerminationCondition.optimal else 'infeasible') +
                  ' L: ' + str(rround(self.local_LB)) + ' U: ' + str(rround(self.global_UB)) + ' Base: ' + str(round(pe.value(self.obj.expr))))
            print('MIP proposes a configuration of:')
            self.print_active_units()
            return True
        else:
            print('MIP infeasible')
            return False

    def solve_set_cover_MIP(self, tee=False):
        # create a mapping of weights on those binaries
        weights = {}
        weights_covered = {cid: 1 for cid in self.covered_units}
        weights_not_covered = {
            cid: len(self.covered_units) + 1 for cid in self.not_covered_units}
        weights.update(weights_covered)
        weights.update(weights_not_covered)
        # use as objective for MILP
        for obj in self.component_objects(ctype=pe.Objective, active=True):
            obj.deactivate()
        self.del_component('init_LOA_obj')
        self.init_LOA_obj = pe.Objective(expr=sum(cid.find_component_on(
            self) * weights[cid] for cid in self.switchable_units), sense=pe.maximize)
        for o in itervalues(self.units):
            o.apply_MIP()
            o.deactivate_oa_cuts()
        results = self.pyomo_results = self.solve(using='local_MIP', tee=tee)
        for o in itervalues(self.units):
            # reactivate OA cuts, otherwise adding another cut will fail
            o.activate_oa_cuts()
        self.solver_time.mip += get_time(results)
        if results.solver.termination_condition is TerminationCondition.optimal:
            print('LOA initialization MIP proposes a configuration of:')
            self.print_active_units()
            return True
        else:
            print('LOA initialization MIP infeasible')
            return False

    def do_bound_contraction(self, **kwargs):
        """Attempts to shrink variable bounds

        Solves a series of bounding MIPs (could also be changed to do LPs)
        in order to determine better bounds for McCormick envelope variables.

        Kwargs:
            verbosity (int): how verbose the output should be.
                0 => suppress all except initial start notification

                1 => notify when a variable bound is contracted

                2 => 1 + notify when no improvement found

                3 => 2 + notify when problem infeasible
        """
        verbosity = kwargs.pop('verbosity', 0)
        print('Performing bound contraction')
        for o in itervalues(self.units):
            getattr(o, 'apply_MIP', doNothing)()

        def varlist():
            """Fetch variables that need to be bound-contracted.
            For now, I will manually determine the variables
            """
            for unit in itervalues(self.units):
                vgen = getattr(unit, 'get_vars_to_bound', doNothing)()
                if vgen is not None:
                    for tup in vgen:
                        # if tup is not a tuple, insert a None for the
                        # block_bounds
                        if isinstance(tup, tuple):
                            yield tup
                        else:
                            yield tup, None

        # Set up upper bounding problem constraints
        self._activate_standard_objective()
        if self.global_UB < float('inf'):
            try:
                self.del_component('global_UB_constr')
            except:
                pass
            self.global_UB_constr = pe.Constraint(
                expr=self.obj.expr <= self.global_UB)

        count = 0
        # Do the bound contraction for each variable
        isInfeasible = False
        # improved bounds: lower bounds, lower block bounds, upper bounds,
        # upper block bounds
        improved_bounds = (0, 0, 0, 0)
        for var, block_bounds in varlist():
            if isInfeasible:
                break
            if var.is_indexed():
                for indx in var._index:
                    if var[indx].is_fixed():
                        continue
                    count += 1
                    result = self._do_bound_contraction(
                        var[indx], block_bounds, verbosity)
                    if result:
                        improved_bounds = tuple(
                            map(sum, zip(improved_bounds, result)))
                    else:
                        isInfeasible = True
                        break
            else:
                if var.is_fixed():
                    continue
                count += 1
                result = self._do_bound_contraction(
                    var, block_bounds, verbosity)
                if result:
                    improved_bounds = tuple(
                        map(sum, zip(improved_bounds, result)))
                else:
                    isInfeasible = True
                    break

        if not isInfeasible:
            print('Ran bound contraction on {} variables'.format(count))
            print('Improved {} LB, {} block LB, {} UB, and {} block UB'.format(
                *improved_bounds))
        else:
            self._global_LB = float('inf')
            return False

        # Remove upper bounding constraint
        try:
            self.del_component('global_UB_constr')
        except:
            pass

        # reconstruct the linear envelopes
        self._reconstruct_envelopes()
        # refresh block bound definitions
        for o in itervalues(self.units):
            o._enforce_block_bounds()
        return True

    def _reconstruct_envelopes(self):
        """Reconstruct the linear envelopes

        Iterates through the list of units and calls the
        'reconstruct_envelopes' function for each unit, if it is defined.
        If subclasses wish to affect which envelopes are reconstrcuted,
        they can override this function.

        TODO only reconstruct the affected linear envelopes

        Returns:
            None
        """
        for o in itervalues(self.units):
            getattr(o, 'reconstruct_envelopes', doNothing)()

    def _do_bound_contraction(self, vardata, block_bounds, verbosity=1):
        """Performs bound contraction on a single _VarData

        Args:
            vardata (_VarData): single variable to bound contract
            block_bounds (dict): dictionary containing active disjunction
                bounds.
                - key: vardata.local_name
                - value: dictionary with
                -- key: bound type, 'ub' or 'lb'
                -- value: bound value or None
            verbosity (int, optional): see do_bound_contraction

        Returns:
            None
        """
        # deactivate other objectives
        for obj in itervalues(self.component_map(
                ctype=pe.Objective, active=True)):
            obj.deactivate()
        lbresult = self._bound_contract('lb', vardata, block_bounds, verbosity)
        if lbresult:
            lb_improved, lb_bb_improved = lbresult
        else:
            return None
        ubresult = self._bound_contract('ub', vardata, block_bounds, verbosity)
        if ubresult:
            ub_improved, ub_bb_improved = ubresult
        else:
            return None
        return lb_improved, lb_bb_improved, ub_improved, ub_bb_improved

    def _bound_contract(self, bnd_type, vardata, block_bounds, verbosity):
        """Contracts a single bound for a single _VarData

        Args:
            bnd_type (str): upper or lower bound, either 'ub' or 'lb'
            vardata (_VarData): single variable to bound contract
            block_bounds (dict): see _do_bound_contraction
            verbosity (int, optional): see do_bound_contraction

        Returns:
            False if infeasible; otherwise, tuple of the number of bounds and
            block bounds improved during the call

        Raises:
            ValueError: if equip_exists variable cannot be found corresponding
            to a _VarData
        """
        # Delete the old bounding objective, if it exists
        try:
            self.del_component('_bounding_obj')
        except:
            pass

        # set the sense and sign for the bounding objective
        if bnd_type == 'lb':
            sense = pe.minimize
            sgn = 1
        elif bnd_type == 'ub':
            sense = pe.maximize
            sgn = -1
        else:
            raise ValueError('Unknown bound of type ' + bnd_type)

        bound_improved = 0  # number of bounds improved
        block_bound_improved = 0  # number of block bounds improved

        # set the bounding objective and solve
        self._bounding_obj = pe.Objective(expr=vardata, sense=sense)
        results = self.pyomo_results = self.solve(using='bounding')

        # record solution time
        self.solver_time.bounding += get_time(results)

        # process bounding problem result
        if results.solver.termination_condition is\
                TerminationCondition.optimal:
            # bounding problem returned a feasible bound
            # check to see if it is tighter than the current value
            new_bnd = pe.value(self._bounding_obj.expr)
            if new_bnd * sgn > getattr(vardata, bnd_type) * sgn:
                if verbosity > 0:
                    print(
                        '{} tightened for {} from {} to {}'.format(
                            bnd_type.upper(),
                            vardata.name,
                            getattr(vardata, bnd_type),
                            new_bnd))
                # vardata.setlb(new_bnd) or .setub(new_bnd)
                getattr(vardata, 'set' + bnd_type)(new_bnd)
                bound_improved += 1
            else:
                if verbosity > 1:
                    print('No improvement found for ' +
                          vardata.name)
        elif results.solver.termination_condition is \
                TerminationCondition.infeasible:
            # bounding problem was infeasible. MIP is infeasible.
            print('Bounding problem for {} of {} was infeasible. Problem '
                  'is infeasible.'.format(
                      bnd_type.upper(), vardata.name))
            return False
        else:
            print('Warning: Failed to solve {} problem for {}'.format(
                bnd_type.upper(), vardata.name))

        # If block_bounds is not None, then calculate the active unit bound
        if block_bounds is not None:
            # get the equip_exists associated with the vardata
            var_parent = vardata.parent_block()
            exists = getattr(var_parent, 'equip_exists', None)
            while exists is None:
                var_parent = var_parent.parent_block()
                if var_parent is None:
                    raise ValueError(
                        'Could not find equip_exists variable among parents of'
                        ' {}'.format(vardata.name))
                exists = getattr(var_parent, 'equip_exists', None)
            exists.fix(1)
            results = self.pyomo_results = self.solve(using='bounding')
            exists.unfix()
            self.solver_time.bounding += get_time(results)
            if results.solver.termination_condition is\
                    TerminationCondition.optimal:
                # check to see if there's a tighter bound possible
                new_bnd = pe.value(self._bounding_obj.expr)
                block_bounds.setdefault(vardata.local_name,
                                        {'lb': None, 'ub': None})
                old_bnd = block_bounds[vardata.local_name][bnd_type]
                if old_bnd is None or new_bnd * sgn > old_bnd * sgn:
                    if verbosity > 0:
                        print(
                            '{} block_bound tightened for {} from {} to {}'.
                            format(
                                bnd_type.upper(),
                                vardata.name,
                                old_bnd,
                                new_bnd))
                    # update bound value
                    block_bounds[vardata.local_name][bnd_type] = new_bnd
                    block_bound_improved += 1
                else:
                    if verbosity > 1:
                        print('No improvement found for ' +
                              vardata.name)
            elif results.solver.termination_condition is \
                    TerminationCondition.infeasible:
                # bounding problem was infeasible. MIP is infeasible.
                print('Block bounding problem for {} of {} was infeasible. '
                      'Unit is fathomed.'.format(
                          bnd_type.upper(), vardata.name))
                # TODO this means that the unit can be fathomed
                var_parent.fathom()
            else:
                print('Warning: Failed to solve {} problem for {}'.format(
                    bnd_type.upper(), vardata.name))

        return bound_improved, block_bound_improved

    def solve_global_NLP(self):
        """Solves the nonlinear program for a fixed flowsheet configuration

        This NLP is used as part of the GLOA algorithm to assess a given
        flowsheet configuration.

        The solver in this case is a global optimization solver, where rigorous
        bounds are guaranteed.

        If the global NLP converges, the flowsheet configuration is feasible.
        An integer cut is added to eliminate the configuration from further
        consideration, since its optimum has been found. If the objective
        value is better than the best incumbent feasible solution, then the
        upper bound is updated with the new value.

        If the global NLP is unable to converge within the given time limit,
        then as much information as possible is retrieved from the results:
         1) if a feasible solution was found, it is evaluated for consideration
         as an improved upper bound. A modified integer cut is also added to
         reflect the lower bound found by the global solver for the given
         configuration.

         2) if no feasible solution was found but an improved lower bound
         was found for the problem, then a modified integer cut is applied to
         reflect the new lower bound.

         3) if no feasible solution or improved lower bound are found, then
         no changes are applied to the problem, but the modeler is given a
         warning that more time may need to be devoted to the global NLP
         solver.

        If the global NLP finds the given configuration provably infeasible,
        then an integer cut is added to eliminate the configuration from
        further consideration.

        Returns:
            None

        Raises:
            NotImplementedError: if an unexpected solver solution status is
                found
        """
        self._activate_standard_objective()
        self._apply_UB_constr()
        for o in itervalues(self.units):
            getattr(o, 'apply_NLP', doNothing)()
        results = self.pyomo_results = self.solve(
            using='global_NLP', load_solutions=False)
        self.solver_time.global_nlp += results.solver.time
        print('Global NLP took ' + str(round(results.solver.time, 2)) + 's')
        if results.solver.termination_condition == TerminationCondition.optimal:
            self.solutions.load_from(results)
            # update upper bound
            self._global_UB = min([self.global_UB, pe.value(self.obj.expr)])
            self.add_integer_cut()
            print('NLP subproblem converged to result of ' +
                  str(pe.value(self.obj.expr)))
        elif results.solver.termination_condition == TerminationCondition.maxTimeLimit:
            # time limit exceeded. If solution exists, load it and update upper
            # bound. Also update lower bound for this problem.
            if results.solution.status == SolutionStatus.feasible:
                # Found a feasible solution. Update upper bound. Add
                # lower-bound integer cut.
                self.solutions.load_from(results)
                if pe.value(self.obj.expr) < self.global_UB:
                    self.solutions.store_to(self.pyomo_results)
                    self._best_found_result = self.pyomo_results
                self._global_UB = min(
                    [self.global_UB, pe.value(self.obj.expr)])
                self.add_mod_integer_cut(results.problem.lower_bound)
                print('NLP subproblem found feasible result of ' +
                      str(pe.value(self.obj.expr)) + ' with LB: ' + str(results.problem.lower_bound))
            elif results.solution.status == SolutionStatus.unknown:
                # Do not know solution status.
                if results.problem.lower_bound > self.global_LB:
                    # update lower bound
                    print('NLP subproblem gives lower bound of: ' +
                          str(results.problem.lower_bound))
                    self.add_mod_integer_cut(results.problem.lower_bound)
                elif results.problem.lower_bound > float('-inf'):
                    print('Warning: no feasible solution or improved lower bound found for NLP subproblem. You may need to allocate more time to the Global NLP solver.')
                else:
                    print('Warning: no feasible solution or lower bound found for NLP subproblem. You may need to allocate more time to the Global NLP solver.')
            else:
                raise NotImplementedError(
                    'Unexpected solver solution status: ' + results.solution.status)
        elif results.solver.termination_condition == TerminationCondition.infeasible:
            self.add_integer_cut()
            print('NLP subproblem proven infeasible.')

        print('NLP Subproblem: ' + (str(round(pe.value(self.obj.expr))) if results.solver.termination_condition ==
                                    TerminationCondition.optimal else 'incomplete') + ' L: ' + str(round(self.global_LB)) + ' U: ' + str(round(self.global_UB)))

    def add_integer_cut(self, tmp=False):
        """Adds an integer cut to the model, eliminating a particular set of
        binary variable values from the space of feasible combinations.

        Args:
            tmp (bool, optional): whether the integer cut should be considered
                'temporary', such as when using LOA to initialize GLOA
        """
        # Add to different ConstraintList depending on whether this is a
        # temporary cut.
        int_cuts = 'int_cuts' if not tmp else 'tmp_int_cuts'
        # Get the appropriate ConstraintList
        getattr(self, int_cuts)\
            .add(  # Add the new integer cut constraint
                # I use the tuple notation here because otherwise I sometimes
                # have a non-fixed bound error from Pyomo
                (1,  # sum(1-y for y^k==1) + sum(y for y^k==0) >= 1
                 sum(1 - o.equip_exists
                     for o in itervalues(self.units)
                     if hasattr(o, 'equip_exists') and
                     abs(o.equip_exists.value - 1.0) <= 1E-3
                     ) +
                 sum(o.equip_exists
                     for o in itervalues(self.units)
                     if hasattr(o, 'equip_exists') and
                     abs(o.equip_exists.value) <= 1E-3
                     ),
                 None
                 ))

    def add_mod_integer_cut(self, sublb, tmp=False):
        """Adds a modified integer cut to the model

        Args:
            sublb (float): the lower bound to associate with the current
                process configuration
            tmp (bool, optional): whether the integer cut should be considered
                'temporary', such as when using LOA to initialize GLOA
        """
        int_cuts = 'int_cuts' if not tmp else 'tmp_int_cuts'
        # Get the appropriate ConstraintList
        getattr(self, int_cuts).add((
            0,
            self.obj.expr - (
                (sublb - self._global_LB) *
                (1 -
                 sum(o.equip_exists
                     for o in itervalues(self.units)
                     if hasattr(o, 'equip_exists') and
                     abs(o.equip_exists.value) <= 1E-3
                     ) -
                 sum(1 - o.equip_exists
                     for o in itervalues(self.units)
                     if hasattr(o, 'equip_exists') and
                     abs(o.equip_exists.value - 1.0) <= 1E-3
                     )
                 ) +
                self._global_LB
            ),
            None
        ))

    def clear_tmp_int_cuts(self):
        """Clears the ConstraintList holding temporary integer cuts

        Returns:
            None
        """
        self.tmp_int_cuts.clear()

    def do_self_proj_cut_gen(self):
        """Performs the self-projection cut generation

        Details for this algorithm are in Francisco's 2016 paper.
        DOI: 10.1016/j.compchemeng.2016.04.017

        Returns:
            None
        """
        for o in itervalues(self.units):
            exists = getattr(o, 'equip_exists', None)
            if exists is not None:
                exists.fix()

        relevant_units = [o for o in self.all_units() if hasattr(
            o, 'generate_cut_gen_problem')]
        # relevant_units = [o for o in itervalues(self.units) if isinstance(o, (Splitter))]
        probs = map(lambda x: x.generate_cut_gen_problem(), relevant_units)

        def solve_cut_gen_problem(cg_prob):
            if cg_prob is None:
                return None, None
            results = self.solve(using='self_proj_cut_gen', model=cg_prob)
            return cg_prob, results
        results = map(solve_cut_gen_problem, probs)

        """ TODO: enable when Pyro-BARON fixed
        Below is a valid implementation of parallel solution with Pyro,
        and it actually works up until the point where it fails because
        the temporary file naming system gums up. So, when that gets fixed,
        I can use this code again.
        """
        # mgr = SolverManagerFactory('pyro')

        # def queue_probs(cg_prob):
        #   if cg_prob is None:
        #       return None, None
        #   globalsolver = SolverFactory('baron')
        #   globalsolver.options['CplexLibName'] = "/opt/ibm/ILOG/CPLEX_Studio1263/cplex/bin/x86-64_linux/libcplex1263.so"
        #   # globalsolver.options['MaxTime'] = 180
        #   globalsolver.options['MaxTime'] = 30
        #   globalsolver.options['allowipopt'] = 0
        #   globalsolver.options['EpsA'] = 0.01
        #   globalsolver.options['EpsR'] = 0.0001
        #   handler = mgr.queue(cg_prob, opt=globalsolver, tee=True)
        #   return cg_prob, handler
        # queued_probs = map(queue_probs, probs)
        # mgr.wait_all()

        # def fetch_results(prob):
        #   cg_prob, handler = prob
        #   if cg_prob is None:
        #       return None, None
        #   else:
        #       return cg_prob, mgr.get_results(handler)
        # results = map(fetch_results, queued_probs)

        def apply_cuts(unit, cg_prob, results):
            if cg_prob is None:
                return
            self.solver_time.self_proj_cut += results.solver.time
            print('Cut gen for ' + unit.local_name + ' took ' +
                  str(round(results.solver.time, 2)) + 's')
            if results.solver.termination_condition is TerminationCondition.optimal:
                unit.apply_self_proj_cut(cg_prob)
            else:
                print(
                    'Unable to solve cut generation problem in time for ' + unit.local_name)

        map(apply_cuts, relevant_units, *zip(*results))

    def add_oa_cut(self):
        """Adds an outer approximation cut at the current solution point

        Advances the internal counter for OA iterations, then calls
        the 'add_oa_cut' function on each of the flowsheet units where it is
        defined.

        Note that an equality relaxation with penalized slack variables
        implementation is used in order to accommodate nonlinear forms.

        Returns:
            None
        """
        self.oa_iter += 1
        for unit in itervalues(self.units):
            unit.add_oa_cut(self.oa_iter)

    def deactivate_oa_cuts(self):
        """Deactivates the outer approximation cuts in the model

        Returns:
            None
        """
        for unit in itervalues(self.units):
            getattr(unit, 'deactivate_oa_cuts', doNothing)()

    def print_active_units(self):
        """Prints a list of the active units in the flowsheet

        Returns:
            None
        """
        unit_tuples = (
            (o.unit_name, o.equip_exists.value) for o in
            itervalues(self.units) if hasattr(o, 'equip_exists') and
            abs(pe.value(o.equip_exists) - 1.0) <= 1E-6
        )
        for name, val in sorted(unit_tuples):
            print(name,)

    def print_all_units(self):
        """Prints a list of all flowsheet units with a binary indicator

        0/1 indicator of whether it is active

        Returns:
            None
        """
        unit_tuples = (
            (o.unit_name, o.equip_exists.value) for o in
            itervalues(self.units) if hasattr(o, 'equip_exists')
        )
        for name, val in sorted(unit_tuples):
            from math import copysign
            print(name, copysign(round(val), 1))

    def register_logic_constraint(self, **kwargs):
        """Registers a logical constraint with the superstructure

        Kwargs:
            one_of (iterable): Registers a logic constraint where one of the
            UnitModel objects in the iterable must be active

        Returns:
            None

        Raises:
            NotImplementedError: the logical key is not recognized
        """
        for key in iterkeys(self.logic):
            term = kwargs.pop(key, ())
            if term:  # not empty
                self.logic[key].add(term)
        for key in iterkeys(kwargs):
            raise NotImplementedError('Unknown logical key: ' + key)

    def build_logic(self):
        """Builds the ConstraintList for superstructure logical constraints

        Scans through items stored in self.logic and adds constraints for
        each of the entries.

        Returns:
            None
        """
        logic_block = getattr(self, 'logic_con', None)
        if logic_block is None:
            logic_block = self.logic_con = pe.ConstraintList()
        for term in self.logic.one_of:
            logic_block.add(
                expr=sum(unit.equip_exists for unit in term) == 1)
        for term in self.logic.at_least_one_of:
            logic_block.add(expr=sum(unit.equip_exists for unit in term) >= 1)
        for term in self.logic.all_of:
            if not isinstance(term, tuple):
                term = (term,)
            for unit in term:
                logic_block.add(expr=unit.equip_exists >= 1)

    def display_costs(self):
        for o in itervalues(self.units):
            cost = getattr(o, 'equip_cost', None)
            exists = getattr(o, 'equip_exists', Bunch(value=1))
            if cost is not None:
                if abs(exists.value - 1.0) <= 1E-6:
                    print(o.unit_name, cost.value)

    def _get_slack_vars(self):
        for o in itervalues(self.units):
            for v in o.get_slack_variables():
                yield v

    def display_slacks(self):
        def compare(x, y):
            if x.value > y.value:
                return 1
            elif x.value < y.value:
                return -1
            else:
                return 0
        print('Active slack penalties')
        for v in sorted(filter(lambda x: x.value > 0, self._get_slack_vars()), cmp=compare):
            print(v.name, v.value)

    def reload_best_found(self):
        if self._best_found_result is None:
            print('Warning: no feasible solution had been found.')
            return
        self.solutions.load_from(self._best_found_result,
                                 ignore_invalid_labels=True)
