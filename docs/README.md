# IDAES Models Documentation #

Last updated: 2016-12-06, Dan Gunter <dkgunter@lbl.gov>

## Overview

This subdirectory is for documentation.

It uses [Sphinx](http://sphinx-doc.org) as the framework.

## Building ##

To build HTML pages:

    make html

## Viewing ##

You can open the webpage locally with your browser. The root 
index.html file is under `_build/html/`. So, e.g., on a Mac you
would type:

    open _build/html/index.html 

You can follow the names of the pages from the table of contents
in each file, starting with 'index.rst',
to get to the appropriate documentation page (ending in
.rst) for a given module. For example: `index.rst` -> `models.rst` -> `superstructure.rst`, where you 
would find/add documentation for the `superstructure_synthesis` module.

There is no magical or automatic association of file names to Python modules and packages, but of course
it is helpful to make it roughly follow the package's structure.

## Jupyter Notebooks ##

Example Jupyter notebooks are stored under "examples/notebooks". 
To include them in the documentation, you can build either restructured text (".rst") files or
full HTML representations. To create these representations, use the Jupyter
"nbconvert" tool. For example:

    jupyter nbconvert --to html examples/notebooks/HelloWorldExample.ipynb

See the file `examples.rst` for the Sphinx include directives. Because Sphinx is able 
to [cross-reference](http://www.sphinx-doc.org/en/1.5/markup/inline.html) between documents,
it may be simpler to put all examples in this file and simply refer to them from within
other parts of the documentation.
