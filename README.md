[![CircleCI](https://circleci.com/gh/IDAES/models.svg?style=shield&circle-token=784380e2298994c6f56411e58d76c3867ec85c0b)](https://circleci.com/gh/IDAES/models)
[![codecov](https://codecov.io/gh/IDAES/models/branch/master/graph/badge.svg?token=vwsYNm2Rdv)](https://codecov.io/gh/IDAES/models)

# Sample solvent model

The purpose of this repository is to get a standard skeleton for documenting
(also running, testing) models.

## Installation


Then, you can run `setup.py` directly or use `pip`. In either case, run the commands from the top-level directory of the repository (directory with this file).

* Using "setup.py":

        python setup.py install
        
    * to install locally, so changes are immediately picked up:
        
            python setup.py develop

* Using "pip":

        pip install .   # note the "." at the end
        
    * to install locally, so changes are immediately picked up:
    
            pip install -e .    # note the "." at the end
            
## Installing solvers
 
Ipopt is a open-source solver that is often used. Prebuilt binaries for Ipopt are available from the link [here](http://ampl.com/products/solvers/open-source/#ipopt).

   

## Testing

The unittests and other tests are under the `tests` directory. To run,

         nose2
         
Note: the `nose2`  command, by default, runs all  functions starting woth `test` in all the Python files starting in `test` that it can find below the current directory. It is part of the [nose2](https://nose2.readthedocs.org/en/latest/) framework.

To get slightly more detailed output, add a verbose flag:

        nose2 -v
    

## Running a notebook

There are example Jupyter notebook(s) in the `examples/` directory. To run them, you should invoke Jupyter on a Notebook file (these end in the extension `.ipynb`, which stands for "IPython Notebook").

    jupyter notebook examples/run-mea-model.ipynb

This should start up a notebook server and then pop up a tab or window in your default web browser showing the Notebook. For more information on how to use Jupyter, see the "Help" menu in the Notebook window itself, and the extensive documentation on the [Jupyter website](https://jupyter.org).

## Documentation

Documentation is in the `docs` directory. It is generated with Sphinx, which uses a Makefile. To build documentation, use the "docs" subcommand of the setup.py script:

    python setup.py docs

Alternatively you can run the Makefile yourself with:

    cd docs
    make html

The main page will be in `_build/html/index.html`.

## License

See file `LICENSE` in this repository.
