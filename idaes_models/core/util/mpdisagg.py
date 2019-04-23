# coding=utf-8
"""
Utility functions for implementing
multi-parametric disaggregation lower bounds

Implementation of lower-bounding formulation given in:
Kolodziej, S., Castro, P. M., & Grossmann, I. E. (2013). Global
    optimization of bilinear programs with a multiparametric
    disaggregation technique. Journal of Global Optimization,
    57(4), 1039–1063. http://doi.org/10.1007/s10898-012-0022-1
"""
from __future__ import division
from pyomo.environ import Constraint, RangeSet, Reals, Binary, NonNegativeReals, Var, Param, Set
from .mccormick import _concat as cc

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def setup_multiparametric_disagg(b, name, minPow, maxPow, *sets):
    setattr(b, name, Constraint(*sets))
    digits = b.digits = RangeSet(0, 9)
    powers = b.powers = RangeSet(minPow, maxPow)
    b.minPow = Param(initialize=minPow)
    b.maxPow = Param(initialize=maxPow)
    skl = sets + (digits, powers)
    sl = sets + (powers,)
    b.y_hat = Var(*skl, domain=Reals)
    b.w_slack = Var(*sets, domain=Reals)
    b.dig_active = Var(*skl, domain=Binary)
    b.x_slack = Var(*sets, domain=NonNegativeReals)
    b.eq_x = Constraint(*sets)
    b.eq_y = Constraint(*sl)
    b.eq_y_hat_lb = Constraint(*skl)
    b.eq_y_hat_ub = Constraint(*skl)
    b.eq_z_sum = Constraint(*sl)
    mc = b.w_slack_mc = Set(initialize=['lb1', 'lb2', 'ub1', 'ub2'])
    sets_mc = sets + (mc,)
    b.eq_w_slack = Constraint(*sets_mc)

def add_mpDisagg_cut(b, name, indx, w, x, y, exists):
    constr = getattr(b, name)
    z = b.dig_active
    b.x_slack[indx].setub(10 ** b.minPow)
    constr.add(indx, expr=w == sum(10 ** l * k * b.y_hat[indx, k, l] for (k, l) in b.digits * b.powers) + b.w_slack[indx])
    b.eq_x.add(indx, expr=x == sum(10 ** l * k * z[indx, k, l] for (k, l) in b.digits * b.powers) + b.x_slack[indx])
    for l in b.powers:
        b.eq_y.add(cc(indx, l), expr=y == sum(b.y_hat[indx, k, l] for k in b.digits))
        b.eq_z_sum.add(cc(indx, l), expr=sum(z[indx, k, l] for k in b.digits) == 1)
        for k in b.digits:
            b.eq_y_hat_lb.add(cc(indx, k) + (l,), expr=y.lb * z[indx, k, l] <= b.y_hat[indx, k, l])
            b.eq_y_hat_ub.add(cc(indx, k) + (l,), expr=b.y_hat[indx, k, l] <= y.ub * z[indx, k, l])
    for bnd in b.w_slack_mc:
        if bnd == 'lb1':
            b.eq_w_slack.add(cc(indx, bnd), expr=b.w_slack[indx] >= y.lb * b.x_slack[indx])
        elif bnd == 'ub1':
            b.eq_w_slack.add(cc(indx, bnd), expr=b.w_slack[indx] <= y.ub * b.x_slack[indx])
        elif bnd == 'lb2':
            b.eq_w_slack.add(cc(indx, bnd), expr=b.w_slack[indx] >= (y - y.ub) * 10 ** b.minPow + y.ub * b.x_slack[indx])
        elif bnd == 'ub2':
            b.eq_w_slack.add(cc(indx, bnd), expr=b.w_slack[indx] <= (y - y.lb) * 10 ** b.minPow + y.lb * b.x_slack[indx])
