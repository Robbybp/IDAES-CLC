'''
IDAES framework

This is skeleton code that will take in models and do something

'''
import types
import argparse
__author__ = "You-Wei Cheah"
__version__ = "1.0.0"


def set_entry_point(func_name):
    """Parameters
    func_name : str
       This is the name of the entry point in a model
    """
    assert isinstance(func_name, types.FunctionType)
    entry_point = func_name


def run(args):
    """Parameters
    args : argparse Namespace
       This is the arguments to be used by the entry point.
    """
    assert isinstance(args, argparse.Namespace)
    entry_point_args = args
