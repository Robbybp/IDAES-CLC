"""
This module provides utility functions and classes related to Pyomo variables
"""
from __future__ import division
from six.moves import range
from six import iteritems
import pyomo.environ as pe

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

class SliceVar(object):
    """
    This class provides a way to pass sliced variables
    into other utility functions
    """
    def __init__(self, var, slices):
        # raise DeprecationWarning('The SliceVar class is deprecated.')
        self.var = var
        self.slices = slices
        # Slices is a dictionary with the format {index: value}
        # So {2: 5, 4: 6} will mean that a call to the SliceVar x with
        # x[1, 3] will return var[1, 5, 3, 6]

    def __getitem__(self, key):
        if key is None:
            key = tuple()
        for pos, value in sorted(iteritems(self.slices)):
            key = key[:pos - 1] + (value,) + key[pos - 1:]
        if len(key) == 0:
            return self.var[None]
        else:
            return self.var[key]

    def __getattr__(self, attr):
        getattr(self.var, attr)

    @property
    def _index(self):
        raise NotImplementedError('I need to overwrite the return value of this to provide the slices of the index set that apply to this SliceVar')

class Wrapper(object):
    """
    This class provides a wrapper around a Pyomo Component so that a Block does
    not try to attach and construct it.

    TODO: this might be a good place for a weakref implementation
    """
    def __init__(self, obj):
        self._obj = obj

    def __getattr__(self, attr):
        return getattr(self._obj, attr)

    def __eq__(self, other):
        """
        Override == operator to test for equlity of the wrapped object instead of the
        wrapper.
        """
        return self._obj == other

def wrap_var(obj):
    """
    Provides a Wrapper around the obj if it is not a constant value; otherwise,
    returns the constant value.
    """
    if pe.is_constant(obj):
        return pe.value(obj)
    else:
        return Wrapper(obj)

def unwrap_var(obj):
    """
    Unwraps the wrapper, if one exists.
    """
    if isinstance(obj, Wrapper):
        return obj._obj
    else:
        return obj

def none_if_empty(tup):
    """Returns None if passed an empty tuple
    This is helpful since a SimpleVar is actually an IndexedVar with
    a single index of None rather than the more intuitive empty tuple.
    """
    if tup is ():
        return None
    else:
        return tup

def lb(expr, block_bounds={}):
    """Returns the lower bound of the expression or variable """
    return _get_bound(expr, 'lb', block_bounds)

def ub(expr, block_bounds={}):
    """Returns the upper bound of the expression or variable """
    return _get_bound(expr, 'ub', block_bounds)

def _get_bound(expr, bound_type, block_bounds={}):
    """Returns the bound of the expression or variable """
    # --- TODO: This is a hack to enable support for floats or ints
    from pyomo.core.base.numvalue import is_constant as const, value
    if const(expr):
        return value(expr)
    # ---
    # if expr.is_constant():
    #     return expr.value
    from pyomo.core.base.var import _GeneralVarData
    if isinstance(expr, _GeneralVarData):
        bound_entry = block_bounds.get(expr.local_name, {'lb': None, 'ub': None})
        block_bound = bound_entry[bound_type]
        if block_bound is not None:
            if bound_type == 'lb':
                return max(block_bound, getattr(expr, bound_type))
            elif bound_type == 'ub':
                return min(block_bound, getattr(expr, bound_type))
        return getattr(expr, bound_type)
    from pyomo.core.base.expr import _SumExpression
    if isinstance(expr, _SumExpression):
        bnd = sum(
            _get_bound(expr._args[i], bound_type, block_bounds) *
            expr._coef[i]
            for i in range(len(expr._args)) if expr._coef[i] > 0) + \
            sum(
                _get_bound(
                    expr._args[i], _invert_bound(bound_type), block_bounds) *
                expr._coef[i]
                for i in range(len(expr._args)) if expr._coef[i] < 0) + \
            expr._const
        return bnd
    from pyomo.core.base.expr_coopr3 import _ProductExpression
    if isinstance(expr, _ProductExpression):
        if len(expr._numerator) == 1 and len(expr._denominator) == 0:
            return expr._coef * _get_bound(expr._numerator[0], bound_type if expr._coef >= 0 else _invert_bound(bound_type), block_bounds)
        # else: (We don't support the expression)
    # else: (We don't recognize the expression)
    raise NotImplementedError(
        'Cannot determine {} for unrecognized expression {}'
        .format(bound_type, expr))

def _invert_bound(bound_type):
    if bound_type == 'lb':
        return 'ub'
    elif bound_type == 'ub':
        return 'lb'
    else:
        raise ValueError('Unknown bound type ' + bound_type)

def is_fixed_by_bounds(expr, block_bounds={}):
    if abs(ub(expr, block_bounds) - lb(expr, block_bounds)) <= 1E-14:
        return True
    else:
        return False

def tighten_var_bound(var, newbound):
    """Tightens the variable bounds for one variable

    This function does not blindly apply the bounds passed in. Rather, it
    evaluates whether the new proposed bounds are better than the existing
    variable bounds.

    Args:
        var (_VarData): single Pyomo variable object (not indexed like Var)
        newbound (tuple): a tuple of (new lower bound, new upper bound)

    Returns:
        None
    """
    lb, ub = var.bounds
    newlb, newub = newbound
    if lb is None or (newlb is not None and newlb > lb):
        var.setlb(newlb)
    if ub is None or (newub is not None and newub < ub):
        var.setub(newub)

def tighten_block_bound(vardata, newbound, block_bounds):
    bb = block_bounds.get(vardata.local_name, {'lb': None, 'ub': None})
    newlb, newub = newbound
    if bb['lb'] is None or (newlb is not None and newlb > bb['lb']):
        bb['lb'] = newlb
    if bb['ub'] is None or (newlb is not None and newub < bb['ub']):
        bb['lb'] = newub

def min_lb(expr1, expr2):
    return min(
        lb(expr1) * lb(expr2),
        lb(expr1) * ub(expr2),
        ub(expr1) * lb(expr2),
        ub(expr1) * ub(expr2))

def max_ub(expr1, expr2):
    return max(
        lb(expr1) * lb(expr2),
        lb(expr1) * ub(expr2),
        ub(expr1) * lb(expr2),
        ub(expr1) * ub(expr2))

def min_lbb(expr1, expr2, block_bounds):
    return min(
        lb(expr1, block_bounds) * lb(expr2, block_bounds),
        lb(expr1, block_bounds) * ub(expr2, block_bounds),
        ub(expr1, block_bounds) * lb(expr2, block_bounds),
        ub(expr1, block_bounds) * ub(expr2, block_bounds))

def max_ubb(expr1, expr2, block_bounds):
    return max(
        lb(expr1, block_bounds) * lb(expr2, block_bounds),
        lb(expr1, block_bounds) * ub(expr2, block_bounds),
        ub(expr1, block_bounds) * lb(expr2, block_bounds),
        ub(expr1, block_bounds) * ub(expr2, block_bounds))

def tighten_mc_var(vardata, x_expr, y_expr, block_bounds):
    tighten_var_bound(
        vardata, (min_lb(x_expr, y_expr), max_ub(x_expr, y_expr)))
    if block_bounds:
        tighten_block_bound(vardata,
                            (min_lbb(x_expr, y_expr, block_bounds),
                             max_ubb(x_expr, y_expr, block_bounds)),
                            block_bounds)
