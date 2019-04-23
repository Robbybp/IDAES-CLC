"""
This class provides a family of wrapper functions that
apply a different nonlinear function depending on whether
a Pyomo or CasADi object is being passed to it.
"""
try:
    import casadi as cas
except:
    pass
    # Does nothing if import not successful, but calling these functions will
    # fail
import pyomo.environ as pyo

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def exp(x):
    if isinstance(x, cas.SX):
        return cas.exp(x)
    else:
        return pyo.exp(x)

def log(x):
    if isinstance(x, cas.SX):
        return cas.log(x)
    else:
        return pyo.log(x)
