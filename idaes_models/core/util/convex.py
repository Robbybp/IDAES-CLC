# coding=utf-8
"""
Utility functions for generating envelopes for convex nonlinear functions
"""
from __future__ import division
from .var import ub, lb, is_fixed_by_bounds, \
    tighten_var_bound, tighten_block_bound
from pyomo.environ import Var, NonNegativeReals, Constraint, RangeSet, Set, Binary
from functools import partial
from .mccormick import squish_concat

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def try_eval(f, x):
    try:
        return f(x)
    except:
        return None

def _setup_convex(b, nsegs, indx):
    if not hasattr(b, 'sets'):
        if indx is None:
            b.sets = Set()
        else:
            b.sets = Set(dimen=len(indx))
    if not hasattr(b, 'overest'):
        if indx is not None:
            b.overest = Constraint(b.sets)
        else:
            b.overest = Constraint()
    if not hasattr(b, 'segs'):
        b.segs = RangeSet(nsegs)
        b.segs_m1 = RangeSet(nsegs - 1, within=b.segs)
    else:
        if nsegs != len(b.segs):
            raise ValueError('We do not currently support different segment counts within the same set of McCormick relaxations.')
    if not hasattr(b, 'x_defn'):
        if indx is not None:
            b.x_defn = Constraint(b.sets, ['lb', 'ub'])
        else:
            b.x_defn = Constraint(['lb', 'ub'])
    if not hasattr(b, 'delta'):
        if indx is not None:
            b.delta = Var(b.sets, b.segs, domain=NonNegativeReals, initialize=0)
        else:
            b.delta = Var(b.segs, domain=NonNegativeReals, initialize=0)
    if not hasattr(b, 'delta_defn'):
        if indx is not None:
            b.delta_defn = Constraint(b.sets, b.segs, ['lb', 'ub'])
        else:
            b.delta_defn = Constraint(b.segs, ['lb', 'ub'])
    if not hasattr(b, 'omega'):
        if indx is not None:
            b.omega = Var(b.sets, b.segs_m1, domain=Binary)
        else:
            b.omega = Var(b.segs_m1, domain=Binary)
    if not hasattr(b, 'omega_exist'):
        if indx is not None:
            b.omega_exist = Constraint(b.sets)
        else:
            b.omega_exist = Constraint()
    if not hasattr(b, 'w'):
        # define bilinear w = y * x
        if indx is not None:
            b.w = Var(b.sets, domain=NonNegativeReals)
        else:
            b.w = Var(domain=NonNegativeReals)
    if not hasattr(b, 'underest'):
        if indx is not None:
            b.underest = Constraint(b.sets, ['lb', 'ub'])
        else:
            b.underest = Constraint(['lb', 'ub'])
    if not hasattr(b, 'w_defn'):
        if indx is not None:
            b.w_defn = Constraint(b.sets, ['lb1', 'ub1', 'lb2', 'ub2'])
        else:
            b.w_defn = Constraint(['lb1', 'ub1', 'lb2', 'ub2'])

def add_convex_relaxation(b, z, x, f_expr, df_expr, nsegs, indx, exists,
                          block_bounds={}, bound_contract=None):
    """Constructs a linear relaxation to bound a convex equality function

    Args:
        b (Block): PyOMO block in which to generate variables and constraints
        z (Expression): PyOMO expression for the convex function output
        x (Expression): PyOMO expression for the convex function input
        f_expr (function): convex function
        df_expr (function): function giving first derivative of convex function
            with respect to x
        exists (_VarData): Variable corresponding to existence of unit
        block_bounds (dict, optional): dictionary describing disjunctive
            bounds present for variables associated with current block

    Returns:
        None
    """
    lbb = partial(lb, block_bounds=block_bounds)  # lower block bound
    ubb = partial(ub, block_bounds=block_bounds)  # upper block bound

    # Define constants
    f_lb = try_eval(f_expr, lb(x))  # z value at lb of x
    f_lbb = try_eval(f_expr, lbb(x))  # z value at lbb of x
    f_ub = try_eval(f_expr, ub(x))  # z value at ub of x
    f_ubb = try_eval(f_expr, ubb(x))  # z value at ubb of x
    df_lb = try_eval(df_expr, lb(x))  # dz/dx value at lb of x
    df_lbb = try_eval(df_expr, lbb(x))  # dz/dx value at lbb of x
    df_ub = try_eval(df_expr, ub(x))  # dz/dx value at ub of x
    df_ubb = try_eval(df_expr, ubb(x))  # dz/dx at ubb of x

    # Tighten bounds on z
    if bound_contract == 'monotonic_increase':
        tighten_var_bound(z, (f_lb, f_ub))
        tighten_block_bound(z, (f_lbb, f_ubb), block_bounds)
    elif bound_contract == 'monotonic_decrease':
        tighten_var_bound(z, (f_ub, f_lb))
        tighten_block_bound(z, (f_ubb, f_lbb), block_bounds)
    else:
        # TODO do nothing for now, but could also solve optimization problem
        # to determine bounds
        pass

    _setup_convex(b, nsegs, indx)

    y = exists  # Alias equipment existence binary as 'y'

    if not is_fixed_by_bounds(x, block_bounds=block_bounds):
        a = (ubb(x) - lbb(x)) / nsegs  # segment length
        b.x_defn.add(squish_concat(indx, 'lb'), expr=x >= lb(x) + sum(b.delta[squish_concat(indx, s)] for s in b.segs) + (lbb(x) - lb(x)) * y)
        b.x_defn.add(squish_concat(indx, 'ub'), expr=x <= ub(x) + sum(b.delta[squish_concat(indx, s)] for s in b.segs) + (lbb(x) - ub(x)) * y)
        b.overest.add(indx, expr=z <= f_lbb + sum((f_expr(lbb(x) + a * s) - f_expr(lbb(x) + a * (s - 1))) / a * b.delta[squish_concat(indx, s)] for s in b.segs) + (ub(z) - f_lbb) * (1 - y))
        for s in b.segs:
            b.delta[squish_concat(indx, s)].setub(a)
            if s < nsegs:
                b.delta_defn.add(squish_concat(indx, s, 'lb'), expr=a * b.omega[squish_concat(indx, s)] <= b.delta[squish_concat(indx, s)])
            if s > 1:
                b.delta_defn.add(squish_concat(indx, s, 'ub'), expr=b.delta[squish_concat(indx, s)] <= a * b.omega[squish_concat(indx, s - 1)])
            else:
                b.delta_defn.add(squish_concat(indx, s, 'ub'), expr=b.delta[squish_concat(indx, s)] <= a * y)
        if nsegs > 1:
            b.omega_exist.add(indx, expr=b.omega[squish_concat(indx, 1)] <= y)

        tighten_var_bound(b.w[indx], (lb(x), ub(x)))

        b.underest.add(squish_concat(indx, 'lb'), expr=z >= df_lb * (x - (lb(x) + (lbb(x) - lb(x)) * y)) + (df_lbb - df_lb) * (b.w - lbb(x) * y) + f_lb + (f_lbb - f_lb) * y)
        b.underest.add(squish_concat(indx, 'ub'), expr=z >= df_ub * (x - (ub(x) + (ubb(x) - ub(x)) * y)) + (df_ubb - df_ub) * (b.w - ubb(x) * y) + f_ub + (f_ubb - f_ub) * y)

        # Add inequalities for overestimators at the lb and ub of x,
        # and Glover envelope for w = y * x
        b.w_defn.add(squish_concat(indx, 'lb1'), expr=b.w >= y * lbb(x))
        b.w_defn.add(squish_concat(indx, 'ub1'), expr=b.w <= y * ubb(x))
        b.w_defn.add(squish_concat(indx, 'lb2'), expr=b.w <= x - (1 - y) * lb(x))
        b.w_defn.add(squish_concat(indx, 'ub2'), expr=b.w >= x - (1 - y) * ub(x))
    else:
        # x is effectively fixed due to its bounds. Impose z = f(x).
        b.underest.add(squish_concat(indx, 'lb'), expr=z >= f_lbb + (lb(z) - f_lbb) * (1 - y))
        b.overest.add(indx, expr=z <= f_ubb + (ub(z) - f_ubb) * (1 - y))
        b.x_defn.add(squish_concat(indx, 'lb'), expr=x >= lb(x) + (lbb(x) - lb(x)) * y)
        b.x_defn.add(squish_concat(indx, 'ub'), expr=x <= ub(x) + (ubb(x) - ub(x)) * y)
        # TODO: use z.fix(f_lb) instead?
