"""
This file sets the Python path so tests can be run without installing the package.
In test files, e.g. `test_foo.py` in this directory, add the following
import statement before other imports of modules in this package:

        from .context import setup_path; setup_path()
"""
import os, sys

def setup_path():
    sys.path.insert(0, os.path.abspath('..'))
    sys.path.insert(0, os.path.abspath('../src/autogenerate'))
