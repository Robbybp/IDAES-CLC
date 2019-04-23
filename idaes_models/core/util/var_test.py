from __future__ import division
from pyomo.environ import value

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def assert_var_equal(test_case, var, expected_val, tolerance):
    test_case.assertIs(
        abs(value(var) - expected_val) <= tolerance, True,
        msg="Value {} not within {} of expected value {}".format(
            value(var),
            tolerance,
            expected_val))

def value_correct(var, expected_val, tolerance):
    return abs(value(var) - expected_val) <= tolerance
