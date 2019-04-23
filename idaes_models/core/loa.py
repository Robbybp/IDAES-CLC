"""Provides functions to perform
local and global logic-based outer approximation
"""
from __future__ import division
import pyomo.environ as pe
from six import itervalues


__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def do_LOA(m, tol=1E-15, iterlim=50, custom_init=None, add_min_flows=True):
    """
    Performs Logic-based Outer Approximation on the FlowsheetModel. The user
    can specify a custom initialization routine to use by passing in a function
    that accepts a FlowsheetModel.

    Args:
        m (FlowsheetModel): flowsheet model
        tol (float, optional): absolute tolerance on bound gap
        iterlim (int, optional): iteration limit for LOA
        custom_init (function, optional): function that accepts a flowsheet
            model as an argument and initializes LOA

    Returns:
        None
    """
    # Initialize LOA
    if custom_init is not None:
        custom_init(m)
    else:
        init_LOA(m, add_min_flows=add_min_flows)

    # Start iterating
    loa_iter = 1
    while m.local_LB + tol < m.global_UB and loa_iter <= iterlim:
        if loa_iter > 1:
            print('')
        print('LOA Iteration {}'.format(loa_iter))
        # return m
        if not m.solve_local_MIP():
            print('LOA terminates on MIP infeasibility')
            break  # if MIP is infeasible, stop.
        last_best_feasible = m.global_UB
        if not m.solve_local_NLP(add_min_flows=add_min_flows):
            print('LOA terminates on NLP infeasibility')
            break  # if NLP is infeasible, stop.
        nlp_improved = m.global_UB < last_best_feasible
        if not nlp_improved:
            print('LOA terminates on NLP non-improvement.')
            break  # if NLP did not find an improvement, stop.
        m.add_oa_cut()
        m.add_integer_cut(tmp=True)
        loa_iter += 1

    if m.local_LB + tol >= m.global_UB and m.global_UB < float('inf'):
        print('LOA algorithm found feasible configuration with objective value {} and non-rigorous LB: {}'.format(round(m.global_UB), round(m.local_LB)))
    else:
        print('LOA algorithm did not converge within iteration limit. LB: ' + str(round(m.local_LB)) + ' UB: ' + str(round(m.global_UB)) + ' after {} iterations'.format(loa_iter - 1))
    print('')

def init_LOA(m, iterlim=5, add_min_flows=True):
    """Performs initialization of the LOA algorithm by solving set covering NLP
    subproblems so that OA linearizations exist for all the major process
    units.

    This work is based upon prototyping done by Eloy Fernandez at CMU.

    Args:
        m (FlowsheetModel): flowsheet
        iterlim (int, optional): iteration limit for initialization

    Returns:
        Boolean: True if set covering successful. False otherwise.
    """
    m.switchable_units = frozenset(pe.ComponentUID(unit.equip_exists) for unit in itervalues(m.units) if hasattr(unit, 'equip_exists'))
    m.covered_units = set()
    m.not_covered_units = set(m.switchable_units)
    iter_count = 1
    while m.not_covered_units and iter_count <= iterlim:
        # print('Covered: ', m.covered_units)
        # print('Not covered: ', m.not_covered_units)
        if not m.solve_set_cover_MIP():
            print('Set covering MIP infeasible. Overall problem infeasible.')
            return False
        # m.solvers.local_NLP.outlev = 3  # use for troubleshooting CONOPT
        if m.solve_local_NLP(add_min_flows=add_min_flows):
            print('Solved initialization NLP. Adding OA cuts.')
            m.add_oa_cut()
            active_units = frozenset(pe.ComponentUID(unit.equip_exists) for unit in itervalues(m.units) if hasattr(unit, 'equip_exists') and abs(pe.value(unit.equip_exists) - 1.0) <= 1E-6)
            m.covered_units.update(active_units)
            m.not_covered_units.difference_update(active_units)
        else:
            # TODO could be infeasible due to bad NLP subproblem
            # initialization. Need to make this more robust
            print('Unable to solve NLP.')
        m.add_integer_cut(tmp=True)
        iter_count += 1
    if m.not_covered_units:
        # Iteration limit was hit without a full covering of all major units.
        print('Iteration limit reached for set covering initialization.')
        return False
    return True

def do_GLOA(m, tol=1E-15, iterlim=1, contract_bounds=False, do_self_proj=False):
    """Performs Global Logic-based Outer Approximation on the FlowsheeetModel

    Args:
        m (FlowsheetModel): flowsheet model
        tol (float, optional): absolute tolerance on bound gap
        iterlim (int, optional): iteration limit for GLOA
        contract_bounds (bool, optional): flag to activate bound contraction
        do_self_proj (bool, optional): flag to activate self projection cuts

    Returns:
        None
    """
    m.clear_tmp_int_cuts()
    m.deactivate_oa_cuts()
    gloa_iter = 1
    m.use_UB_constr = True
    while m.global_LB + tol < m.global_UB and gloa_iter <= iterlim:
        if gloa_iter > 1:
            print('')  # new line
        print('GLOA Iteration ' + str(gloa_iter))
        if contract_bounds:
            if not m.do_bound_contraction():
                print('GLOA terminates on bounding infeasibility.')
                break
        if not m.solve_MIP():
            print('GLOA terminates on MIP infeasibility')
            break
        if do_self_proj:
            raise NotImplementedError('Self projection cuts are currently not supported, until separation problem generation is fixed.')
            m.do_self_proj_cut_gen()
        m.solve_global_NLP()
        gloa_iter += 1
    print('')

    if m.global_LB + tol >= m.global_UB and m.global_UB < float('inf'):
        print('GLOA algorithm found optimal configuration with objective value {} and LB: {}'.format(round(m.global_UB), round(m.global_LB)))
    else:
        print('GLOA algorithm did not converge within iteration limit. LB: ' + str(round(m.global_LB)) + ' UB: ' + str(round(m.global_UB)) + ' after {} iterations'.format(gloa_iter - 1))
