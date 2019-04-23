"""
A simple GUI viewer for Pyomo models.
"""
from __future__ import division
from __future__ import print_function
__author__ = "John Eslick"
__version__ = "1.0.0"

import warnings
try:
    from PyQt5.QtWidgets import QApplication
except ImportError:
    try:
        from PyQt4.QtGui import QApplication
    except:
        warnings.warn("Cannot import PyQt")

from qt_viewer.variable_browser import *
import datetime

class QtModelViewer():

    def __init__(self, model=None, qt_args=[]):
        self.model = model
        self.app = QApplication(qt_args)
        print("Created viewer")
        print(datetime.datetime.time(datetime.datetime.now()))

    def show_vars(self):
        self.w = VariableViewer(self.model)
        self.app.exec_()


if __name__ == "__main__":
    """
    Test viewer
    """
    from pyomo.environ import *

    model = ConcreteModel()
    model.T = Var(initialize=300, domain=PositiveReals,
                  bounds = (250,500), doc="Temperature (K)")
    model.P = Var(initialize=101325, domain=PositiveReals,
                  bounds = (1e5,1e6), doc="Pressure (Pa)")
    v = QtModelViewer(model)
    v.show_vars()
    print(value(model.T))
    print(value(model.P))
