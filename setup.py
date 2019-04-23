import os, sys
import subprocess
from setuptools import setup, find_packages, Extension, Command

sphinx_requires = ['Sphinx', 'alabaster', 'sphinxcontrib-napoleon',
'sphinx-rtd-theme']
other_requires = [
    'ipython',
    'nose',
    'nose2',
    'pandas',
    'six',
    'snowballstemmer',
    'networkx',
    'pyomo>=4.3',
    'bunch'
]
class BuildDocsCommand(Command):
    """Convenience command for building the documentation.
    """
    description = 'Build Sphinx documentation'
    user_options = [('target=', 't', 'Target for "make" in docs directory')]

    VALID_TARGETS = ('html', 'dirhtml', 'singlehtml', 'pickle', 'json', 'htmlhelp',
                     'qthelp' ,'devhelp', 'epub', 'latex', 'latexpdf', 'text', 'man',
                     'changes', 'linkcheck', 'doctest', 'coverage', 'gettext')

    def initialize_options(self):
        self.target = 'html'

    def finalize_options(self):
        if self.target not in self.VALID_TARGETS:
            print('ERROR: Invalid documentation target "{}". Valid options: {}'.format(
                self.target, ', '.join(self.VALID_TARGETS)))
            sys.exit(1)

    def run(self):
        try:
            os.chdir('docs')
            subprocess.call(['make', self.target])
        finally:
            os.chdir('..')

class BuildPropCommand(Command):
    """Convenience command for building the physical properties C library.
    """
    description = 'Build C physical properties library, used for ExternalFunctions in tests.'
    user_options = [('adol=', 'a', 'Path to ADOL-C install, e.g: /usr/local/src/pyadolc/PACKAGES/ADOL-C/ADOL-C'),
    ('asl=', 'a', 'Path to ASL install, e.g: /usr/local/src/Ipopt-3.12.5/ThirdParty/ASL')]

    def initialize_options(self):
        # Initializing to Docker paths.
        self.adol = '/usr/local/src/pyadolc/PACKAGES/ADOL-C/ADOL-C'
        self.asl = '/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL'

    def finalize_options(self):
        pass

    def run(self):
        try:
            os.chdir('idaes_models/properties/physical/MEA_Simple')
            
            # Take out previous builds:
            command = ['make','clean']            
            program = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            program.wait()
            
            # Build: 
            command = ['make','ADOL={0}'.format(self.adol),'ASL_BUILD={0}'.format(self.asl)]            
            program = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            program.wait()

            if os.path.exists('phys_prop.so'):
                print ("Successful physical properties build.")
            else:
                print ("Physical properties build encountered errors.")
                # Output result of build:
                if program is not None:
                    print(program.communicate()[0])    
        finally:
            os.chdir('../../../..')

setup(
    name="idaes_models",
    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    version="0.1.0",
    install_requires=['pyomo'] + sphinx_requires + other_requires,
    author="IDAES Project",
    author_email="idaes@host.org",
    maintainer="IDAES Project",
    url="https://github.com/idaes/models",
    license="MIT",
    description="IDAES project models",
    long_description="Simulation models for the IDAES project",
    keywords=["IDAES"],
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    entry_points = {
        'console_scripts': [
            'wrap = idaes_models.python_c_wrapper.autogenerate.wrap:main',
        ],
    },
    scripts=['scripts/run_model'],
    cmdclass={
        'docs': BuildDocsCommand,
        'phys_prop': BuildPropCommand
    }
)
