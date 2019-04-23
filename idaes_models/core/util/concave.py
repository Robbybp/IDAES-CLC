# coding=utf-8
"""
Utility functions for implementing
piecewise linear underestimators for concave univariate expressions.

Implementation of the delta formulation from:
Bergamini, M. L., Grossmann, I., Scenna, N., & Aguirre, P. (2008).
    An improved piecewise outer-approximation algorithm for the global
    optimization of MINLP models involving concave and bilinear terms.
    Computers and Chemical Engineering, 32, 477–493.
    http://doi.org/10.1016/j.compchemeng.2007.03.011
"""
from __future__ import division
from pyomo.environ import RangeSet, Param, NonNegativeReals, Var, Constraint, Binary, Set
from .mccormick import _concat as cc, squish_concat
from .var import ub, lb, is_fixed_by_bounds, tighten_var_bound, \
    tighten_block_bound
from functools import partial

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def add_concave_linear_underest(b, name, nsegs, x, f, f_expr, *sets, **kwargs):
    exists = kwargs.get('exists', 1.0)
    bb = kwargs.get('block_bounds', {})
    if not isinstance(sets, tuple):
        sets = (sets,)
    segs = b.segs = RangeSet(nsegs)
    segs_m1 = b.segs_m1 = RangeSet(nsegs - 1)
    sets_segs = sets + (segs,)
    sets_segs_m1 = sets + (segs_m1,)

    def calc_lengths(m, *i):
        i = i if i else None  # change empty tuple to None
        return (ub(x[i], bb) - lb(x[i], bb)) / nsegs
    seg_length = b.seg_length = Param(*sets, initialize=calc_lengths)

    delta = b.delta = Var(*sets_segs, domain=NonNegativeReals)
    for i_seg in delta._index:
        if not isinstance(i_seg, tuple):
            # if we're getting a singleton (not defined over a set)
            i = None
        else:
            i = i_seg[:-1]
        delta[i_seg].setub(seg_length[i])
    w = b.w = Var(*sets_segs_m1, domain=Binary)

    def x_defn(m, *i):
        i = i if i else None  # change empty tuple to None
        return x[i] == lb(x[i], bb) + sum(delta[cc(i, s)] for s in segs) - \
            (1 - exists) * (lb(x[i], bb) - lb(x[i]))
        # TODO: this restricts the upper bound of x[i] when the unit does not exist. This is fine if x[i] = 0 when the unit doesn't exist, but is a hack otherwise.
    b.x_defn = Constraint(*sets, rule=x_defn)

    def f_lb(m, *i):
        if len(i) == 0:
            i = None
        if abs(seg_length[i]) < 1E-10:
            return f[i] >= f_expr(*cc(i, lb(x[i], bb))) - \
                (1 - exists) * (f_expr(*cc(i, lb(x[i], bb))) -
                                f_expr(*cc(i, lb(x[i]))))
        return f[i] >= f_expr(*cc(i, lb(x[i], bb))) + \
            sum((f_expr(*cc(i, lb(x[i], bb) + s * seg_length[i])) -
                 f_expr(*cc(i, lb(x[i], bb) + (s - 1) * seg_length[i]))
                 ) / seg_length[i] * delta[cc(i, s)]
                for s in segs) - \
            (1 - exists) * (f_expr(*cc(i, lb(x[i], bb))) -
                            f_expr(*cc(i, lb(x[i]))))
    b.f_lb = Constraint(*sets, rule=f_lb)

    def delta_lb(m, *i_seg):
        # split out the segment from the sets
        i, seg = i_seg[:-1], i_seg[-1]
        return seg_length[i] * w[i, seg] <= delta[i, seg] if seg < nsegs else Constraint.NoConstraint
    b.delta_lb = Constraint(*sets_segs, rule=delta_lb)

    def delta_ub(m, *i_seg):
        i, seg = i_seg[:-1], i_seg[-1]
        return delta[i, seg] <= seg_length[i] * w[i, seg - 1] if seg > 1 else Constraint.NoConstraint
    b.delta_ub = Constraint(*sets_segs, rule=delta_ub)

def try_eval(f, x):
    try:
        return f(x)
    except:
        return None

def _setup_concave(b, nsegs, indx):
    if not hasattr(b, 'sets'):
        if indx is None:
            b.sets = Set()
        else:
            b.sets = Set(dimen=len(indx))
    if not hasattr(b, 'underest'):
        if indx is not None:
            b.underest = Constraint(b.sets)
        else:
            b.underest = Constraint()
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
    if not hasattr(b, 'overest'):
        if indx is not None:
            b.overest = Constraint(b.sets, ['lb', 'ub'])
        else:
            b.overest = Constraint(['lb', 'ub'])
    if not hasattr(b, 'w_defn'):
        if indx is not None:
            b.w_defn = Constraint(b.sets, ['lb1', 'ub1', 'lb2', 'ub2'])
        else:
            b.w_defn = Constraint(['lb1', 'ub1', 'lb2', 'ub2'])

def add_concave_relaxation(b, z, x, f_expr, df_expr, nsegs, indx, exists, block_bounds={}, bound_contract=None):
    lbb = partial(lb, block_bounds=block_bounds)
    ubb = partial(ub, block_bounds=block_bounds)

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

    _setup_concave(b, nsegs, indx)

    y = exists  # Alias equipment existence binary as 'y'

    if not is_fixed_by_bounds(x, block_bounds=block_bounds):
        a = (ubb(x) - lbb(x)) / nsegs  # segment length
        b.x_defn.add(squish_concat(indx, 'lb'), expr=x >= lb(x) + sum(b.delta[squish_concat(indx, s)] for s in b.segs) + (lbb(x) - lb(x)) * y)
        b.x_defn.add(squish_concat(indx, 'ub'), expr=x <= ub(x) + sum(b.delta[squish_concat(indx, s)] for s in b.segs) + (lbb(x) - ub(x)) * y)
        b.underest.add(indx, expr=z >= f_lbb + sum((f_expr(lbb(x) + a * s) - f_expr(lbb(x) + a * (s - 1))) / a * b.delta[squish_concat(indx, s)] for s in b.segs) + (lb(z) - f_lbb) * (1 - y))
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

        if df_lb is not None:
            b.overest.add(squish_concat(indx, 'lb'), expr=z <= df_lb * (x - (lb(x) + (lbb(x) - lb(x)) * y)) + (df_lbb - df_lb) * (b.w - lbb(x) * y) + f_lb + (f_lbb - f_lb) * y)
        if df_ub is not None:
            b.overest.add(squish_concat(indx, 'ub'), expr=z <= df_ub * (x - (ub(x) + (ubb(x) - ub(x)) * y)) + (df_ubb - df_ub) * (b.w - ubb(x) * y) + f_ub + (f_ubb - f_ub) * y)

        # Add inequalities for overestimators at the lb and ub of x,
        # and Glover envelope for w = y * x
        b.w_defn.add(squish_concat(indx, 'lb1'), expr=b.w >= y * lbb(x))
        b.w_defn.add(squish_concat(indx, 'ub1'), expr=b.w <= y * ubb(x))
        b.w_defn.add(squish_concat(indx, 'lb2'), expr=b.w <= x - (1 - y) * lb(x))
        b.w_defn.add(squish_concat(indx, 'ub2'), expr=b.w >= x - (1 - y) * ub(x))
    else:
        # x is effectively fixed due to its bounds. Impose z = f(x).
        b.underest.add(indx, expr=z >= f_lbb + (lb(z) - f_lbb) * (1 - y))
        b.overest.add(squish_concat(indx, 'lb'), expr=z <= f_ubb + (ub(z) - f_ubb) * (1 - y))
        b.x_defn.add(squish_concat(indx, 'lb'), expr=x >= lb(x) + (lbb(x) - lb(x)) * y)
        b.x_defn.add(squish_concat(indx, 'ub'), expr=x <= ub(x) + (ubb(x) - ub(x)) * y)
        # TODO: use z.fix(f_lb) instead?
