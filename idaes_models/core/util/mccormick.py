# coding=utf-8
"""
Utility functions for implementing
piecewise (or standard) McCormick envelopes

Implementation of the McCormick envelope formulation given in:
Misener, R., Thompson, J. P., & Floudas, C. A. (2011). Apogee:
    Global optimization of standard, generalized, and extended pooling
    problems via linear and logarithmic partitioning schemes.
    Computers and Chemical Engineering, 35(5), 876–892.
    http://doi.org/10.1016/j.compchemeng.2011.01.026
"""
from __future__ import division
from pyomo.environ import Set, Constraint, Param, Var, RangeSet, Binary, NonNegativeReals
from .var import lb, ub, is_fixed_by_bounds
from functools import partial

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def squish_concat(a, *b):
    return squish(_concat(a, *b))

def _concat(a, *b):
    if isinstance(a, tuple):
        return a + b
    elif a is None:
        return b
    else:
        return (a,) + b

def squish(tup):
    """Squishes a singleton tuple ('A',) to 'A'

    If tup is a singleton tuple, return the underlying singleton. Otherwise,
    return the original tuple.

    Args:
        tup (tuple): Tuple to squish
    """
    if len(tup) == 1:
        return tup[0]
    else:
        return tup

def _setup_mccormick(b, nsegs, indx):
    if not hasattr(b, 'mccormick_cuts'):
        b.mccormick_cuts = Set(initialize=['lb1', 'lb2', 'ub1', 'ub2'])
    if not hasattr(b, 'sets'):
        if indx is None:
            b.sets = Set()
        elif isinstance(indx, tuple):
            b.sets = Set(dimen=len(indx))
        else:
            b.sets = Set(dimen=1)
    if indx is not None and indx not in b.sets:
        # Presumably I would need to add indices to this set, but
        # apparently it's complaining when I do, so I'm not.
        # b.sets.add(indx)
        pass
    if not hasattr(b, 'constr'):
        if indx is not None:
            b.constr = Constraint(b.sets, b.mccormick_cuts)
        else:
            b.constr = Constraint(b.mccormick_cuts)
    if not hasattr(b, 'segs'):
        b.segs = RangeSet(nsegs)
    else:
        if nsegs != len(b.segs):
            raise ValueError('We do not currently support different segment counts within the same set of McCormick relaxations.')
    if not hasattr(b, 'xbar'):
        if indx is not None:
            b.xbar = Var(b.sets, initialize=0)
        else:
            b.xbar = Var(initialize=0)
    if not hasattr(b, 'delta_y'):
        if indx is not None:
            b.delta_y = Var(b.sets, b.segs, domain=NonNegativeReals, initialize=0)
        else:
            b.delta_y = Var(b.segs, domain=NonNegativeReals, initialize=0)
    if not hasattr(b, 'lbda'):
        if indx is not None:
            b.lbda = Var(b.sets, b.segs, domain=Binary)
        else:
            b.lbda = Var(b.segs, domain=Binary)
    if not hasattr(b, 'x_defn'):
        if indx is not None:
            b.x_defn = Constraint(b.sets, ['lb', 'ub'])
        else:
            b.x_defn = Constraint(['lb', 'ub'])
    if not hasattr(b, 'xbar_defn'):
        if indx is not None:
            b.xbar_defn = Constraint(b.sets, ['lb1', 'lb2', 'ub1', 'ub2'])
        else:
            b.xbar_defn = Constraint(['lb1', 'lb2', 'ub1', 'ub2'])
    if not hasattr(b, 'y_defn'):
        if indx is not None:
            b.y_defn = Constraint(b.sets, ['lb', 'ub'])
        else:
            b.y_defn = Constraint(['lb', 'ub'])
    if not hasattr(b, 'dy_ub'):
        if indx is not None:
            b.dy_ub = Constraint(b.sets, b.segs)
        else:
            b.dy_ub = Constraint(b.segs)
    if not hasattr(b, 'lbda_defn'):
        if indx is not None:
            b.lbda_defn = Constraint(b.sets)
        else:
            b.lbda_defn = Constraint()
    b.xbar._index.add(indx)
    for s in b.segs:
        # b.delta_y._index.add((indx, s))
        b.delta_y.add(squish(_concat(indx, s)))
        # b.lbda._index.add(indx)
        b.lbda.add(squish(_concat(indx, s)))

def add_mccormick_relaxation(b, z, x, y, nsegs, indx, exists, block_bounds={}):
    """Adds McCormick envelopes for a bilinear term z = x * y

    Args:
        b (Block): PyOMO block in which to put constraints and helper variables
        z (Expression): PyOMO expression for the bilinear product
        x (Expression): expression for the bilinear operand to be divided into
            segments for the piecewise case
        y (Expression): expression for the other bilinear operand
        nsegs (integer): number of piecewise segments (normal is 1)
        indx (tuple or singleton): index for an indexed bilinear term
        exists (Var): variable corresponding to equipment existence
        block_bounds (dict, optional): dictionary describing disjunctive
            bounds present for variables associated with current block

    Returns:
        None
    """
    lbb = partial(lb, block_bounds=block_bounds)
    ubb = partial(ub, block_bounds=block_bounds)
    _setup_mccormick(b, nsegs, indx)
    delta_y = b.delta_y
    lbda = b.lbda
    x_defn = b.x_defn
    if not is_fixed_by_bounds(x, block_bounds=block_bounds):
        a = (ubb(x) - lbb(x)) / nsegs  # segment length
        if 'lb1' in b.mccormick_cuts:
            b.constr.add(
                squish(_concat(indx, 'lb1')),
                expr=z >= b.xbar[indx] * lbb(y) + sum((lbb(x) + a * (s - 1)) * delta_y[squish(_concat(indx, s))] for s in b.segs) + (1 - exists) * lb(z))
        if 'lb2' in b.mccormick_cuts:
            b.constr.add(
                squish(_concat(indx, 'lb2')),
                expr=z >= b.xbar[indx] * ubb(y) + sum((lbb(x) + a * s) * (delta_y[squish(_concat(indx, s))] - (ubb(y) - lbb(y)) * lbda[squish(_concat(indx, s))]) for s in b.segs) + (1 - exists) * lb(z))
        if 'ub1' in b.mccormick_cuts:
            b.constr.add(
                squish(_concat(indx, 'ub1')),
                expr=z <= b.xbar[indx] * lbb(y) + sum((lbb(x) + a * s) * delta_y[squish(_concat(indx, s))] for s in b.segs) + (1 - exists) * ub(z))
        if 'ub2' in b.mccormick_cuts:
            b.constr.add(
                squish(_concat(indx, 'ub2')),
                expr=z <= b.xbar[indx] * ubb(y) + sum((lbb(x) + a * (s - 1)) * (delta_y[squish(_concat(indx, s))] - (ubb(y) - lbb(y)) * lbda[squish(_concat(indx, s))]) for s in b.segs) + (1 - exists) * ub(z))
        x_defn.add(squish(_concat(indx, 'lb')), expr=lb(x) + (lbb(x) - lb(x)) * exists + sum(a * (s - 1) * lbda[squish(_concat(indx, s))] for s in b.segs) <= x)
        x_defn.add(squish(_concat(indx, 'ub')), expr=x <= ub(x) + (ubb(x) - ub(x)) * exists + sum(a * s * lbda[squish(_concat(indx, s))] for s in b.segs))
    else:
        # x fixed by bounds
        if 'lb1' in b.mccormick_cuts:
            b.constr.add(
                squish(_concat(indx, 'lb1')),
                expr=z >= b.xbar[indx] * lbb(y) + sum(lbb(x) * delta_y[squish(_concat(indx, s))] for s in b.segs) + (1 - exists) * lb(z))
        if 'ub1' in b.mccormick_cuts:
            b.constr.add(
                squish(_concat(indx, 'ub1')),
                expr=z <= b.xbar[indx] * lbb(y) + sum(ubb(x) * delta_y[squish(_concat(indx, s))] for s in b.segs) + (1 - exists) * ub(z))
    b.xbar_defn.add(squish(_concat(indx, 'lb1')), expr=b.xbar[indx] >= exists * lbb(x))
    b.xbar_defn.add(squish(_concat(indx, 'ub1')), expr=b.xbar[indx] <= exists * ubb(x))
    b.xbar_defn.add(squish(_concat(indx, 'lb2')), expr=b.xbar[indx] <= x - (1 - exists) * lb(x))
    b.xbar_defn.add(squish(_concat(indx, 'ub2')), expr=b.xbar[indx] >= x - (1 - exists) * ub(x))
    b.y_defn.add(squish(_concat(indx, 'lb')), expr=y >= lb(y) + (lbb(y) - lb(y)) * exists + sum(delta_y[squish(_concat(indx, s))] for s in b.segs))
    b.y_defn.add(squish(_concat(indx, 'ub')), expr=y <= ub(y) + (ubb(y) - ub(y)) * exists + sum(delta_y[squish(_concat(indx, s))] for s in b.segs))
    for s in b.segs:
        b.dy_ub.add(squish(_concat(indx, s)), expr=delta_y[squish(_concat(indx, s))] <= (ubb(y) - lbb(y)) * lbda[squish(_concat(indx, s))])
    b.lbda_defn.add(indx, expr=sum(lbda[squish(_concat(indx, s))] for s in b.segs) == exists)

# Depr
def setup_mccormick_cuts(b, name, nsegs, *sets):
    print('Warning: setup_mccormick_cuts is now deprecated in favor of the util.mccormick.add_mccormick_relaxation function, which automatically performs setup.')
    cuts = b.mccormick_cuts = Set(initialize=['lb1', 'lb2', 'ub1', 'ub2'])
    sets_cuts = sets + (cuts,)
    setattr(b, name, Constraint(*sets_cuts))
    if nsegs > 1:
        # If we are doing piecewise McCormick, set up additional constraints and variables.
        b.seg_length = Param(*sets, mutable=True)
        segs = b.segs = RangeSet(nsegs)
        sets_segs = sets + (segs,)
        b.seg_active = Var(*sets_segs, domain=Binary)
        b.delta_y = Var(*sets_segs, domain=NonNegativeReals)
        b.eq_seg_active_sum = Constraint(*sets)
        b.eq_y = Constraint(*sets)
        b.eq_dy_ub = Constraint(*sets)
        b.eq_x_lb = Constraint(*sets)
        b.eq_x_ub = Constraint(*sets)

# Depr
def add_mccormick_cut(b, name, indx, z, x, y, exists, **kwargs):
    print('Warning: add_mccormick_cut is now deprecated in favor of the util.mccormick.add_mccormick_relaxation function, which supports block bounds and equipment existence with piecewise capability.')
    bb = kwargs.get('block_bounds', {})
    constr = getattr(b, name)
    cuts = b.mccormick_cuts
    cc = _concat
    if hasattr(b, 'segs'):
        # Piecewise McCormick implementation
        # TODO enable expressions support and exists support for piecewise
        # TODO enable block bounds for piecewise
        a = b.seg_length[indx] = (x.ub - x.lb) / len(b.segs)
        dy = b.delta_y
        for seg in b.segs:
            dy[indx, seg].setub(y.ub - y.lb)
        seg_active = b.seg_active
        if 'lb1' in cuts:
            constr.add(
                cc(indx, 'lb1'),
                expr=z >= x * y.lb + sum(
                    (x.lb + a * (seg - 1)) * dy[indx, seg]
                    for seg in b.segs))
        if 'lb2' in cuts:
            constr.add(cc(indx, 'lb2'), expr=z >= x * y.ub + sum((x.lb + a * seg) * (dy[indx, seg] - (y.ub - y.lb) * seg_active[indx, seg]) for seg in b.segs))
        if 'ub1' in cuts:
            constr.add(cc(indx, 'ub1'), expr=z <= x * y.lb + sum((x.lb + a * seg) * dy[indx, seg] for seg in b.segs))
        if 'ub2' in cuts:
            constr.add(cc(indx, 'ub2'), expr=z <= x * y.ub + sum((x.lb + a * (seg - 1)) * (dy[indx, seg] - (y.ub - y.lb) * seg_active[indx, seg]) for seg in b.segs))

        b.eq_seg_active_sum.add(indx, expr=sum(seg_active[indx, seg] for seg in b.segs) == 1)
        b.eq_y.add(indx, expr=y == y.lb + sum(dy[indx, seg] for seg in b.segs))
        for seg in b.segs:
            b.eq_dy_ub.add(cc(indx, seg), expr=dy[indx, seg] <= (y.ub - y.lb) * seg_active[indx, seg])
        b.eq_x_lb.add(indx, expr=x.lb + sum(a * (seg - 1) * seg_active[indx, seg] for seg in b.segs) <= x)
        b.eq_x_ub.add(indx, expr=x <= x.lb + sum(a * seg * seg_active[indx, seg] for seg in b.segs))
    else:
        # Non-piecewise McCormick
        if 'lb1' in cuts:
            constr.add(
                squish(cc(indx, 'lb1')),
                expr=z >= lb(x, bb) * y + x * lb(y, bb) -
                lb(x, bb) * lb(y, bb) -
                (max(lb(x, bb) * lb(y, bb), lb(x, bb) * ub(y, bb)) +
                 max(lb(x, bb) * lb(y, bb), ub(x, bb) * lb(y, bb)) -
                 lb(x, bb) * lb(y, bb) - lb(z)) * (1 - exists))
        if 'lb2' in cuts:
            constr.add(
                squish(cc(indx, 'lb2')),
                expr=z >= ub(x, bb) * y + x * ub(y, bb) -
                ub(x, bb) * ub(y, bb) -
                (max(ub(x, bb) * lb(y, bb), ub(x, bb) * ub(y, bb)) +
                 max(lb(x, bb) * ub(y, bb), ub(x, bb) * ub(y, bb)) -
                 ub(x, bb) * ub(y, bb) - lb(z)) * (1 - exists))
        if 'ub1' in cuts:
            constr.add(
                squish(cc(indx, 'ub1')),
                expr=z <= ub(x, bb) * y + x * lb(y, bb) -
                ub(x, bb) * lb(y, bb) -
                (min(ub(x, bb) * lb(y, bb), ub(x, bb) * ub(y, bb)) +
                 min(lb(x, bb) * lb(y, bb), ub(x, bb) * lb(y, bb)) -
                 ub(x, bb) * lb(y, bb) - ub(z)) * (1 - exists))
        if 'ub2' in cuts:
            constr.add(
                squish(cc(indx, 'ub2')),
                expr=z <= x * ub(y, bb) + lb(x, bb) * y -
                lb(x, bb) * ub(y, bb) -
                (min(lb(x, bb) * ub(y, bb), ub(x, bb) * ub(y, bb)) +
                 min(lb(x, bb) * lb(y, bb), lb(x, bb) * ub(y, bb)) -
                 lb(x, bb) * ub(y, bb) - ub(z)) * (1 - exists))
