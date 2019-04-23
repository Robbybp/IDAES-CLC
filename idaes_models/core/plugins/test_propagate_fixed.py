from idaes_models.core.plugins.propagate_fixed import propagate_var_fix
from pyomo.environ import ConcreteModel, Var, Constraint
import unittest
from idaes_models.core.util.misc import category

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

class TestPropagateFixed(unittest.TestCase):
    @category('frequent')
    def test_var_set(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=2)
        m.v3 = Var(initialize=3)
        m.v4 = Var(initialize=4)
        m.c1 = Constraint(expr=m.v1 == m.v2)
        m.c2 = Constraint(expr=m.v2 == m.v3)
        m.c3 = Constraint(expr=m.v3 == m.v4)
        m.v2.fix()
        # check to make sure that all the v's have the same equality set. John
        # had found a logic error.
        propagate_var_fix(m)
        self.assertTrue(m.v1.fixed)
        self.assertTrue(m.v2.fixed)
        self.assertTrue(m.v3.fixed)
        self.assertTrue(m.v4.fixed)
        # m.display()


if __name__ == '__main__':
    unittest.main()
