"""
Base class for property packages
"""
# Chages the divide behavior to not do integer division
from __future__ import division
from __future__ import print_function

# Some more inforation about this module
__author__ = "John Eslick"
__version__ = "1.0.0"

from pyomo.environ import *
from process_base import ProcessBase
import copy
import six
import model_util
import sys
import itertools
import traceback
#from collections import OrderedDict

__all__ = ['PropertyPackage']

class _Prop():
    """
    
    """
    def __init__(self, name, doc="", unit="", bounds=(0,1),
                 domain=PositiveReals, indexSet=None):
        """
        Args:
        name: string name of this property
        func: is the name of function to call to add this prop
        doc: description of this property
        unit: string representing unit of measure
        """
        self.name = name
        self.doc = doc
        self.unit = unit
        self.bounds = bounds
        self.domain = domain
        self.indexSet = indexSet

class _PropCalc():
    """
    This is a class to hold some information about functions that add
    properties to a pyomo block
    
    Attributes:
    ===========
    func: function to do inital value cals and/or set up variables and
          constraints.
    ig_func:
    dep: List of strings for functions that this directly depends on
    provides: List of properties this provides
    """
    def __init__(self, func, ig_func, dep=(), provides=(), constraints=()):
        """
        Args:
        =====
        func: function to add variables/constraints to a property block
        ig_func: function to do initial guess calculations
        dep: is the list of dependencies for the property
        doc: is a string that describes this propery
        unit: is a unit of meassure string
        var_list: is a list of variable names for setting inital guess
        """
        self.func = func
        self.ig_func = ig_func
        self.dep = dep
        self.provides = provides
        if constraints is None:
            constraints = ()
        self.constraints = constraints
        
class _ConstraintModAdd(pyomo.core.base.constraint.IndexedConstraint):
    """
    This makes a slight change to the way the add method works. This
    is meant to be temporary need to check on this
    """
    def add(self, *args, **kwargs):
        """
        The first argument of add should be the index.  If it is a tuple
        of lenght 1 replace by the scalar, otherwise leave unchanged.
        The call the add method of the base class.
        """
        if isinstance(args[0], tuple) and len(args[0]) == 1:
            args = list(args)
            args[0]=args[0][0]
        return pyomo.core.base.constraint.IndexedConstraint\
            .add(self, *args, **kwargs)
        
class PropertyPacakage(ProcessBase):
    """
    Base class for properties package
    """
    def __init__(self, comp, properties=[], T_init=330, P_init=101325, 
                 y_init=None, index_set=((0,),)):
        """
            Also want molar volume as possible initial guess, but will
            set that aside for now.
            
            Args:
            =====
            comp:  A list of strings representing chemical components
            properties: A list of properties to calculate T, P, and y
                        are standard.
            T_init: Initial temperature guess in K
            P_init: Initial pressure guess in Pa
            y_init: Initial moles fraction guess
            indexSet: Optional set to create an index on all properties
                      for discretization mostly.
        """
        ProcessBase.__init__(self, index_set)
        self.comp = comp # all components
        # A set of component pairs, no duplication, where ordering 
        # dosen't matter. (for like binary diffusivites...)
        self.binary_set = []
        for i in self.comp:
            for j in self.comp: 
                if i!=j and (j,i) not in self.binary_set:
                    self.binary_set.append((i,j))
        # Create empty property dictionary and propery set up function
        # list
        self._prop_dict = {}  #Dict of properties/variables available
        self._prop_calc = []  #List of functions to setup property calc
        self._calc_lookup = {}#Dict to lookup function for a property
        # Set the initial guess for standard variables that are always
        # required for a property package.
        self.std_props = ('T', 'P', 'y')
        # Set the list of properties to calculate from user
        # depending on what they specify some may be added later
        self.properties = properties # list of properties to calculate
        # Create the list of standard properties
        self._register_properties()
        # Register the setup functions for the required props (T, P, y)
        self._reg_calc(provides=self.std_props,
                       func=self.standard_vars, 
                       ig_func=self.standard_vars_ig,
                       constraints=None)
        self._reg_initialize(
            'default', func=self._default_init, 
            doc='Fix T, P, and y then solve')
        self.T_init = T_init
        self.P_init = P_init
        self.y_init = y_init
        
    def _default_init(self):
        """
        Default init, simple, should work fine for all property methods
        """
        # copy the block, so can change and solve without touching 
        # original
        blk = self.copy_pyomo()
        # Unfix all the varaiables not sure what may have been fixed in
        # original.  Not worried about active/inactive constraints here
        # should be okay.
        model_util.unfix_all(blk)
        # Fix a set of variables, which should fully specify properties
        model_util.fix_vars(blk, ['T','P','y'])
        # Create a new model that will be used to solve block equations
        tmp_model = ConcreteModel(name='tmp_init')
        # Add block to model
        tmp_model.pb = blk
        # For ipot at least you don't even need objective should have 
        # 0 degrees of freedom, just leave next line commented as 
        # reminder may need it
        #tmp_model.obj = Objective(expr=1)
        # Create the solver
        slv = SolverFactory('ipopt')
        # Solve
        slv.solve(tmp_model)
        # Should check result here and raise exception if didn't solve
        # will need to look into that.
        #
        #
        # Save the results in a dict
        d = self.save_json_static(blk, dict_out=True)
        # Return the dict containting the solved values
        return d
        
    def _block_constraints(self, blk):
        for fi in self._call_list:
            f = self._prop_calc[fi]
            for c in f.constraints:
                if c[1] is not None:
                    indxst = self.index_set + c[1]
                else:
                    indxst = self.index_set
                cnst = _ConstraintModAdd(*indxst, noruleinit=True)
                setattr(blk, c[0], cnst)
                
    def _block_vars(self, blk):
        """
        Create a Pyomo variable for registered properties.  The idea
        is that this will simplify some things like pulling the bounds
        and initial guess from the registered properties, and allow an
        easy way to add an index to the property variables when 
        discretizing a property method.
        """
        for p in self.properties:
            lb = self._prop_dict[p].bounds[0]
            ub = self._prop_dict[p].bounds[1]
            ig = lb + (ub-lb)/4.0
            if self._prop_dict[p].indexSet is not None:
                args = self.index_set + self._prop_dict[p].indexSet
            else:
                args = self.index_set
            kwargs = {
                'initialize':ig,
                'domain':self._prop_dict[p].domain,
                'bounds':(lb,ub),
                'doc':self._prop_dict[p].doc}
            setattr(blk, p, Var(*args, **kwargs))
            
    def _register_properties(self):
        """
        This registers the properties defined in the standard. Specific
        properties packages can also define their own properties, 
        especially intermediate values that are specific to certain
        """
        r = self._reg_prop
        
        r('T', unit='K', domain=PositiveReals, doc='Temperature',
          bounds=(2e2, 1e3))
          
        r('P', unit='Pa', domain=PositiveReals, doc='Pressure', 
          bounds=(1, 1e8))
        
        r('y', unit='mol/mol', domain=NonNegativeReals,
          doc='Mole fraction', bounds=(0,1), indexSet=(self.comp,))
          
        r('V_vap', unit='m^3/mol', domain=PositiveReals, 
          doc='Vapor molar volume', bounds=(1e-4, 1e-1))
        
        r('V_liq_pure', unit='m^3/mol', domain=PositiveReals, 
          doc='Vapor molar volume', bounds=(1e-6, 1e-4), 
          indexSet=(self.comp,))
        
        r('V_liq', unit='m^3/mol',  domain=PositiveReals,
          doc='Liquid molar volume', bounds=(1e-6, 1e-4))
          
        r('gamma', unit='none', domain=PositiveReals,
          doc='Activity coefficient', bounds=(0,1))
        
        r('rho_liq', unit='kg/m^3',  domain=PositiveReals,
          doc='Liquid mass density', bounds=(1e2, 5e3))
        
        r('rho_vap', unit='kg/m^3',  domain=PositiveReals,
          doc='Vapor mass density', bounds=(0,1000))
        
        r('D_liq', unit='m^2/s',  domain=PositiveReals,
          doc="Diffusivity of components in liquid", bounds=(1e-11,1e-8),
          indexSet=(self.comp,))
        
        r('D_vap', unit='m^2/s', domain=PositiveReals,
          doc="Diffusivity of components in vapor", bounds=(1e-7,1e-4),
          indexSet=(self.comp,))
        
        r('D_vap_ij', unit='m^2/s', domain=PositiveReals,
          doc="Binary diffusivity of components in vapor",
          bounds=(1e-7,1e-4), indexSet=(self.binary_set,))           
        
        r('phi', unit='none', domain=PositiveReals,
          doc="Fugacity coefficient", bounds=(0,1))
        
        # Add Gibbs
        
        r('henry', unit='Pa*m^3/mol', domain=PositiveReals,
          doc="Henry's constant", bounds=(100, 1e5),
          indexSet=(self.comp,))
        
        r('cp_mol_vap_ideal_pure', unit='J/mol/K', domain=PositiveReals,
          doc="Constant pressure vapor pure component ideal gas heat "\
              "capacity", bounds=(1,100), indexSet=(self.comp,))
        
        r('cp_mol_vap_ideal', unit='J/mol/K', domain=PositiveReals,
          doc="Constant pressure vapor ideal gass heat capacity",
          bounds=(1,100))
        
        r('cp_mol_vap', unit='J/mol/K', domain=PositiveReals,
          doc="Constant pressur vapor heat capacity", bounds=(1,100))
        
        r('cv_mol_vap_ideal_pure', unit='J/mol/K', domain=PositiveReals,
          doc="Constant volume vapor pure component ideal gas heat "\
              "capacity", bounds=(1,100), indexSet=(self.comp,))
              
        r('cv_mol_vap_ideal', unit='J/mol/K', domain=PositiveReals,
          doc="Constant volume vapor ideal gass heat capacity",
          bounds=(1,100))
          
        r('cv_mol_vap', unit='J/mol/K', domain=PositiveReals, 
          doc="Constant volume vapor heat capacity", bounds=(1,100))
          
        r('cp_mol_liq', unit='J/mol/K', domain=NonNegativeReals,
          doc="Constant pressure liquid heat capacity", bounds=(0,100000))
        
        r('cp_mol_liq_pure', unit='J/mol/K', domain=NonNegativeReals,
          doc="Pure component constant pressure liquid heat capacity", 
          bounds=(0,100000), indexSet=(self.comp,))
          
        r('cv_mol_liq', unit='J/mol/K', domain=NonNegativeReals,
          doc="Constant volume liquid heat capacity", bounds=(0,1000))
          
        # Add Entropy
        # Add Enthalpy
        
        r('dh_vap', unit='J/mol', domain=NonNegativeReals,
          doc="Component heat of vaporization", 
          bounds=(100, 1e6), indexSet=(self.comp,))
        
        r('dh_abs', unit='J/mol', domain=NonNegativeReals,
          doc="Heat of absorption", bounds = (100,1e6), 
          indexSet=(self.comp,))                       
        
        r('MW', unit='kg/mol',  domain=PositiveReals,
          doc='Mixture molecular weight', bounds=(0.001, 0.200))
        
        r('wt_frac', unit='kg/kg', domain=PositiveReals, bounds=(0,1),
          doc='Component weight fractions', indexSet=(self.comp,))
                       
        r('sigma', unit='N/m', domain=PositiveReals,
          doc="Surface tension", bounds=(0,10))
        
        r('sigma_pure', unit='N/m', domain=PositiveReals,
          doc="Pure component surface tension", bounds=(0,10), 
          indexSet=(self.comp,))
        
        r('kc_liq_pure', unit='W/(m*K)', domain=PositiveReals,
          doc="Pure component liquid thermal conductivity", 
          bounds=(1e-1,100), indexSet=(self.comp,))
        
        r('kc_liq', unit='W/(m*K)', domain=PositiveReals,
          doc="Liquid thermal conductivity", bounds=(1e-1,10))
          
        r('kc_vap_pure', unit='W/(m*K)', domain=PositiveReals,
          doc="Pure component vapor thermal conductivity", 
          bounds=(1e-2,1), indexSet=(self.comp,))
          
        r('kc_vap', unit='W/(m*K)', domain=PositiveReals,
          doc="Vapor thermal conductivity", bounds=(1e-2,1)) 
        
        r('P_vap', unit='Pa', domain=NonNegativeReals,
          doc="Vapor pressure of components", bounds=(1,1e5), 
          indexSet=(self.comp,))
        
        r('P_vap_pure', unit='Pa', domain=PositiveReals,
          doc="Vapor pressure of components", bounds=(1,1e5), 
          indexSet=(self.comp,))
          
        r('mu_liq', unit='kg/m/s', domain=PositiveReals,
          doc="Liquid viscosity", bounds=(1e-6, 0.05))
          
        r('mu_vap', unit='kg/m/s', domain=PositiveReals,
          doc="Vapor viscosity", bounds=(1e-7, 1e-4),)
        
        r('mu_vap_pure', unit='kg/m/s', domain=PositiveReals,
          doc="Vapor pure component viscosity", bounds=(1e-7, 1e-4),
          indexSet=(self.comp,))
          
        r('conc_vap_tot', unit='mol/m^3', domain=PositiveReals,
          doc='Vapor concentration/mole density', bounds=(1,1e2))
          
        r('conc_vap', unit='mol/m^3', domain=NonNegativeReals,
          doc='Vapor component concentrations', bounds=(0,1e2), 
          indexSet=(self.comp,))
          
        r('conc_liq_tot', unit='mol/m^3', domain=PositiveReals,
          doc='Liquid concentration/mole density', bounds=(1000,4e5)) 
          
        r('conc_liq', unit='mol/m^3', domain=NonNegativeReals,
          doc='Liquid component concentrations', bounds = (10,1e5), 
          indexSet=(self.comp,))

    def _expand_guess(self, var, var_index=()):
        """
        This function takes an intial guess and if the dimensions don't
        match the index_set epands it out. Either by using a scalar or
        interpolating.
        """
        try: 
            keys = var.keys()
            miss = False
            for i in itertools.product(*(self.index_set+var_index)):
                if i not in keys:
                    miss = True
                    break
            if not miss:
                # All the guesses are there don't need anything
                return var
        except AttributeError:
            # not a dict, so assume it is a scalar in this case assign
            # the same value to all intial guesses and return
            v = {}
            for i in itertools.product(*(self.index_set+var_index)):
                v[i] = var
            return v
        # Now see if an intial guess is indexed but not over index set
        i = keys[0]
        if not isinstance(i, tuple):
            i = (i,)
        if len(i) < len(self.index_set+var_index):
            #assume guess has the var_var index but not index_set
            v = {}
            idxs_len = len(self.index_set)
            for i in itertools.product(*(self.index_set+var_index)):
                i2 = i[idxs_len:]
                if len(i2) == 1:
                    i2 = i2[0]
                v[i] = var[i2]
            return v
        # Now try interpolation, if here var is a dict and has all
        # dimension but some initial guesses are missing
        raise Exception("Interpolatng initial guesses not implimented")
        
    def standard_vars_ig(self, z, blk):
        """
        Create an initial guess for standard T, P, and y variables. All
        other initial guess calcaultion stem from these so need to get
        dimensions right.  Rest should be fine.
        """
        blk.T[z].value = self.T_init[z]
        blk.P[z].value = self.P_init[z]
        for i in self.comp:
            blk.y[z,i].value = self.y_init[z+(i,)]
    
    def standard_vars(self, z, blk):
        """
        No constraints for standard vars.
        """
        pass
        
    def __call__(self, blk, index=0):
        """
        Build Pyomo equation block
        
        Args:
        =====
        blk: Pyomo blk being created
        indx: index of Pyomo block being created, currently no use for it
        """
        self.pyomo_block = blk
        self._block_vars(blk)
        self._block_constraints(blk)
        # Now create constraints
        self.T_init = self._expand_guess(self.T_init)
        self.P_init = self._expand_guess(self.P_init)
        self.y_init = self._expand_guess(
            self.y_init, var_index=(self.comp,))
        for fi in self._call_list:
            f = self._prop_calc[fi]
            for z in itertools.product(*self.index_set):
                f.ig_func(z, blk)
                f.func(z, blk)
        self.pyomo_block = blk
        return blk
    
    def _reg_prop(self, name, unit='', doc='', bounds=(0,1), 
                  domain=PositiveReals, indexSet=None):
        """
        Information about properties that can be calculated by this
          (If there are some properties or variables that are just
           intermediate values don't really need to list them here
           but won't hurt if you do)
        """
        self._prop_dict[name] = _Prop(name, unit=unit, doc=doc, 
                                      bounds=bounds, domain=domain, 
                                      indexSet=indexSet)
    
    def _unreg_prop(self):
        try:
            del(self._prop_dict[name])
        except:
            pass
        
    def _reg_calc(self, provides=(), dep=(), func=None, ig_func=None,
        constraints=()):
        """
        Register information about functions to add property calculations
        """
        self._prop_calc.append(_PropCalc(
            provides=provides, dep=dep, func=func, ig_func=ig_func,
            constraints=constraints))
        f = self._prop_calc[-1]
        i = len(self._prop_calc) - 1
        for p in f.provides:
            self._calc_lookup[p] = i

    def set_props(self, vlist=None):
        """
        Set the property list and figure out the property dependecies 
        and the order in which to call the functions.  The order is 
        important for the initial guess calculation.
        
        Args:
        =====
        vlist: list of properties to calculate, if None use the object's
               vlist, vlist=None is for the dependecy check/ordering
        """
        # If no property list (vlist), use list from object
        if vlist is None:
            vlist = self.properties
        # make sure vlist a set, order doesn't matter and only want 
        # things listed once.
        vlist = set(vlist)
        # Make sure vlist includes required standard vars T, P, and y
        vlist.update(self.std_props)
        # Create a set of functions to call to
        try:
            func_list = set([self._calc_lookup[p] for p in vlist])
        except KeyError, e:
            raise Exception("No function provides {0}".format(e))
            
        # Create a new vlist with all provided properies
        vset = set()
        depset = set()
        for f in func_list:
            vset.update(self._prop_calc[f].provides)
            depset.update(self._prop_calc[f].dep)
        # Expand vlist with dependencies.
        for p in depset:
            vset.add(p)
        # Now need to sort so the dpendencies are done before needed.
        # The order really only affects the initial guesses, although
        # it could also be important if there are varaible shared 
        # between functions that are not registered as properties.
        #
        # Not worried about efficeny, time required is insignificant
        props = set(vset) # vlist is good to go, so make copy to change
        props_added = set() # properties in nl
        prop_order = []
        nl = 1 # angthing not zero is fine
        while props and len(prop_order) != nl:
            nl = len(prop_order) #make sure done get stuck in a loop
            for p in props:
                deps = self._prop_calc[self._calc_lookup[p]].dep
                if props_added.issuperset(deps):
                    props_added.add(p)
                    prop_order.append(p)
            props.difference_update(props_added)
        if props:
            print("Warning you property package has circular"
                  "dependencies. Check the following:")
            print(props)
            # could still work with ig error trapping
            prop_order.extend(props)
        
        self._call_list = []
        for p in prop_order:
            fidx = self._calc_lookup[p]
            if fidx not in self._call_list:
                self._call_list.append(fidx)
        self.properties = vset
         
    def _check_ig(self, l=()):
        """
        Check if inital guess is already available for variable names
        listed in l. Returns a list varialbes not found.
        
        Args:
        l = list of strings that need to be checked for an inital guess
        """
        not_found = []
        for i in l:
            if not hasattr(self.ig, i):
                not_found.append(i)
        return not_found
