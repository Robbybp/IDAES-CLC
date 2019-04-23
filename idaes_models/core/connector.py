from __future__ import division
from pyomo.environ import Connector


# TODO [from Qi]: Is this code used anymore?
def fluid_connector(Fvar=None, Tvar=None, Pvar=None):
    """
    Basic fluid stream connector

    Args:
        Fvar (Var, optional): Pyomo variable corresponding to component flows
        Tvar (Var, optional): Pyomo variable corresponding to temperature
        Pvar (Var, optional): Pyomo variable corresponding to pressure
    """
    c = Connector()
    if Fvar is not None:
        c.add(Fvar, name="F")
    if Tvar is not None:
        c.add(Tvar, name="T")
    if Pvar is not None:
        c.add(Pvar, name="P")
    return c
