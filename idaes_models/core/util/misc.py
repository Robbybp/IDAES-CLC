"""
Miscellaneous utilty functions
"""
import pyutilib.services

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


def doNothing(*args, **kwargs):
    """This function does nothing. It is useful for instances when you want to
    call a function, if it exists. For example: getattr(unit,
    'possibly_defined_function', getNothing)()

    Args:
        *args (anything): accepts any argument
        **kwargs (anything): accepts any keyword arguments

    Returns:
        None
    """
    pass


def get_time(results):
    time = getattr(results.solver, 'time', None)
    if time is not None:
        return time
    time = getattr(results.solver, 'wallclock_time', None)
    if time is not None:
        return time
    time = getattr(results.solver, 'User time', None)
    if time is not None:
        return time
    else:
        raise NotImplementedError('Could not determine solver time')


def rround(n, *args, **kwargs):
    try:
        return round(n, *args, **kwargs)
    except OverflowError:
        return n


def category(*args):
    """Decorator for tests to enable tiered testing.

    Suggested categories:
        1. frequent
        2. nightly
        3. expensive
        4. research

    Args:
        *args (tuple of strings): categories to which the test belongs

    Returns:
        function: Either the original test function or skip
    """
    import os
    import unittest

    if 'idaes_test_level' in os.environ:
        _active_categories = os.environ['idaes_test_level'].strip()
        try:
            # See if the environment variable is a number. If yes, then run all
            # the categories up and including that number.
            test_level = int(_active_categories)
            categories = ('frequent', 'nightly', 'expensive', 'research')
            _active_categories = categories[:test_level]
        except ValueError:
            # For now, only support entry of one category here
            _active_categories = (_active_categories,)
    else:
        _active_categories = ('frequent',)

    if 'idaes_test_exclude' in os.environ:
        _exclude_categories = os.environ['idaes_test_exclude'].strip()
        _exclude_categories = (_exclude_categories,)
    else:
        _exclude_categories = ()

    if any(cat in args for cat in _exclude_categories):
        return unittest.skip('Test categories {} contains one of excluded categories "{}"'.format(sorted(args), sorted(_exclude_categories)))
    elif any(cat in args for cat in _active_categories):
        def wrapper(func):
            return func
        return wrapper
    else:
        return unittest.skip('Test categories {} do not match active categories "{}"'.format(sorted(args), sorted(_active_categories)))

def requires_solver(solver):
    from pyomo.opt import SolverFactory
    import unittest

    if not SolverFactory(solver).available():
        return unittest.skip('Required solver {} is not available.'.format(solver))
    else:
        def wrapper(func):
            return func
        return wrapper


def get_pyomo_tmp_files():
    """
    Make Pyomo write it's temporary files to the current working directory,
    useful for checking nl, sol, and log files for ASL solvers without needing
    to track down the temporary file location.
    """
    pyutilib.services.TempfileManager.tempdir = './'


def hhmmss(sec_in):
    """
    Convert elapsed time in seconds to "d days hh:mm:ss.ss" format.
    This is nice for things that take a long time.
    """
    h = int(sec_in // 3600)
    m = int(sec_in % 3600 // 60)
    s = sec_in % 3600 % 60
    if h < 24:
        hstr = "{0:0>2}".format(h)
    elif h >= 24 and h < 48:
        hstr = "1 day {0:0>2}".format(h % 24)
    else:
        hstr = "{0} days {1:0>2}".format(h // 24, h % 24)
    return "{0}:{1:0>2}:{2:0>5.2f}".format(hstr, m, s)
