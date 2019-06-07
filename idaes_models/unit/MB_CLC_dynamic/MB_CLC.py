#!/,usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:55:31 2019

Moving Bed Reactor Model - based on Kim et al. (2016) and Mondino et al. (2017)
Counter-current (gas phase inlet @ z=0, solid phase inlet @ z=1)
Steady-state

Assumptions: - Ideal Gas behavior
             - Constant solids velocity through the bed
             - Non-isothermal/adiabatic/non-adiabatic/ 
               heat loss through refractory-lined reactor wall
             - Interstitial gas velocity calculated using Ergun Eq. or 
               simplified pressure balance

Change Log (key changes):
    Improved model initialization routine, including initializing gas profile 
        temperatures at solid inlet temperature as thermal mass of solid >>
        thermal mass of gas
    Merged ergun and simplified P eqns as options for pressure drop correlation
    Made updates to the gas energy balance eqn as vg, cp_gas and rho_vap are 
        functions of axial length. Used a thermal flux variable to capture this
    Updated documentation of eqns and variables to improve legibility

Constraint blocks:
    a constraints - Bed geometry (eq_a1 to eq_a2)
    b constraints - Mass balance (eq_b1 to eq_b4)
    c constraints - Flowrate and flux relationships (eq_c1 to eq_c17)
    d constraints - Energy balance (eq_d1 to eq_d10)
    e constraints - Momentum/pressure balance (eq_e1 to eq_e4)
    f constraints - Boundary and outlet conditions (eq_f1 to eq_f10)
    g constraints - Transfer coefficients (eq_g1 to eq_g14)
    h constraints - Initial conditions
    p constraints - Thermophysical properties (eq_p1 to eq_p16)
    r constraints - Reaction kinetic properties (eq_r1 to eq_r6)
"""
from __future__ import division
from __future__ import print_function              
__author__ = "Chinedu Okoli and Anca Ostace"
__version__ = "2.0.0"


from pyomo.environ import Param, Var, Constraint, Expression, Set, Reals, \
                            acos,exp, value, TransformationFactory, \
                            NonNegativeReals
from pyomo.dae import DerivativeVar, ContinuousSet

import time
#import math
#import importlib

# Import IDAES cores
from idaes_models.core import UnitModel, ProcBlock
#import idaes_models.core.util.misc as model_util

from pyomo.opt import SolverFactory

import pdb

#__all__ = ['MB Reactor']

@ProcBlock("MB")
class _MB(UnitModel):
    
    # Lists of all chemical species and reactions included in package
    scomp = ['Fe2O3','Fe3O4','Al2O3']    
    gcomp = ['CO2','H2O','CH4']
    comp = scomp + gcomp
    reaction_idx = [1]
    
    def __init__(self, *args, **kwargs):
        """
        Args:
        dae_method = method to use for calcuating derivatives (default = OCLR)
                    - BFD1 - 1st order backwards finite difference
                    - OCLR - Orthogonal collocation, Lagrange-Radau
                    - OCLL - Orthogonal collocation, Lagrange-Legendre
        fe_set = set of normalised finite element locations
        nfe = number of finite elements for bed discretization (default = 15)
                    (not used if fe_set specified)
        ncp = number of collocation points (OCLR or OCLL only, default = 3)
        press_drop = Pressure drop correlation for superficial velocity calc.
                    - SimplifiedP - simplified pressure correlations
                    - Ergun - Ergun equation

        fe_set_t = set of time finite elements
        nfe_t = number of finite elements in time domain
        ncp_t = number of col points in time domain
        H = time horizon length 


        """     
        
        self.dae_method = kwargs.pop("dae_method","OCLR")

        # why does this create its own fe set instead of
        # using TransformationFactory?
        # if I'm not mistaken, both are used? 
        # why? This seems like a recipe for making a mistake...
        # is this because constraints will be added and removed
        # during initialization?
        self.fe_set = kwargs.pop("fe_set",[])            
        if self.fe_set == []:
            # Get nfe and populate fe_set
            self.nfe = kwargs.pop("nfe", 15)
            for i in range(0,self.nfe+1):
                if i == 0:
                    self.fe_set.append(0.0)
                elif i < self.nfe:
                    self.fe_set.append(float(i/self.nfe))
                else:
                    self.fe_set.append(1.0)
        else:
            # Set nfe from len(fe_set)
            self.nfe = len(self.fe_set) - 1  
            kwargs.pop("nfe", 15)

        # need bounds to construct fe set here
        #self.fe_set_t = kwargs.pop('fe_set_t',[])
        self.horizon = kwargs.pop('horizon',60.0)
        self.nfe_t = kwargs.pop('nfe_t',12)
        self.ncp_t = kwargs.pop('ncp_t',3)
        #if self.fe_set_t == []:
        #    self.nfe_t = kwargs.pop('nfe_t',15)
        #    for i in range(0,self.nfe_t+1):
        #        if i == 0:
        #            self.fe_set_t.append(0.0)
        #        elif i < self.nfe_t:
        #            self.fe_set_t.append(float(i/self.nfe_t))
        #        else: 
        #            self.fe_set_t.append(self.horizon)
        #else:
        #    # if a fe set is provided (which it shouldn't be...)
        #    self.nfe = len(self.fe_set_t) - 1
        #    # what is this line:?
        #    kwargs.pop('nfe',15)

        
            
        print(self.fe_set)  
#        print(self.fe_set_t)
        
        self.ncp = kwargs.pop("ncp",3)

        self.press_drop = kwargs.pop("press_drop","Ergun")
        
        
        UnitModel.__init__(self, *args, **kwargs)

    def build(self, *args, **kwargs):
        self._make_params_props()
        self._make_params()
        self._make_domain()

        self._make_init_values()
        
        self._make_vars()
        self._make_vars_props()
        self._make_derivatives()
        
        self._make_constraints_props()
        self._make_bed_geometry()
        self._mass_balance()
        self._energy_balance()
        self._pressure_balance()
        self._make_bdry_conds()

        self._make_init_conditions()
        
        self._transfer_coefficients()
        

        self._dae_transform()
        self._t_dae_transform()


        
        
    def _make_params_props(self):
        """
        Create the gas and solid properties paramters.
        These would belong to a property package, but have been moved here
        because it somehow improved convergence. Will be moved back out 
        while transitioning to the new framework.
        """        
        # Declare gas and solid component lists
        self.GasList = Set(initialize=self.gcomp)
        self.SolidList = Set(initialize=self.scomp)
#        self.comp = Set(initialize=self.comp)
        self.rxn_idx = Set(initialize=self.reaction_idx)
        
        # Dummy parameter used to scale the reaction rate expression during
        # initialization
        self.scale = Param(within = Reals, default = 1, mutable = True, 
                         doc= 'scale factor for rate rxn')
        self.Tref = Param(default=298.15,
                    doc='Thermodynamic Reference Temperature [K]')
        
        # Mol. weights of solids - units = kg/kmol. ref: NIST webbook
        self.MW = {'CH4':16,'CO2':44,'H2O':18, 
                     'Fe2O3':159.69, 'Fe3O4':231.533, 'Al2O3':101.96}
        
        # Stoichiometric coefficients
        '''Stoichiometric coefficient for each component in each reaction'''
        self.stoic = {'CH4':-1, 'Fe2O3':-12, 'Fe3O4':8, 'CO2':1, 'H2O':2, 'Al2O3':0}

        # Std. heat of formation of components - units = J/(mol comp) - ref: NIST
        self.Hf_comp = {'CH4':-74873.1,'CO2':-393522.4,'H2O':-241826.4,
                        'Fe2O3':-825503.2,'Fe3O4':-1120894,'Al2O3':-1675690}
        # Ideal gas spec. heat capacity parameters(Shomate) of components - ref: NIST webbook
        # Shomate equations from NIST. Parameters A-E are used for cp calcs while A-H are used for enthalpy calc
        # cp_comp = A + B*T + C*T^2 + D*T^3 + E/(T^2) where T = Temperature (K)/1000, and cp_comp = (J/mol.K)
        # H_comp = H - H(298.15) = A*T + B*T^2/2 + C*T^3/3 + D*T^4/4 - E/T + F - H where T = Temp (K)/1000 and H_comp = (kJ/mol)     
        self.cp_param = {
                        ('Al2O3', 'a'): 102.4290,
                        ('Al2O3', 'b'): 38.74980,
                        ('Al2O3', 'c'): -15.91090,
                        ('Al2O3', 'd'): 2.628181,
                        ('Al2O3', 'e'): -3.007551,
                        ('Al2O3', 'f'): -1717.930,
                        ('Al2O3', 'g'): 146.9970,
                        ('Al2O3', 'h'): -1675.690,
                        ('Fe3O4', 'a'): 200.8320000,
                        ('Fe3O4', 'b'): 1.586435e-7,
                        ('Fe3O4', 'c'): -6.661682e-8,
                        ('Fe3O4', 'd'): 9.452452e-9,
                        ('Fe3O4', 'e'): 3.18602e-8,
                        ('Fe3O4', 'f'): -1174.1350000,
                        ('Fe3O4', 'g'): 388.0790000,
                        ('Fe3O4', 'h'): -1120.8940000,
                        ('Fe2O3', 'a'): 110.9362000,
                        ('Fe2O3', 'b'): 32.0471400,
                        ('Fe2O3', 'c'): -9.1923330,
                        ('Fe2O3', 'd'): 0.9015060,
                        ('Fe2O3', 'e'): 5.4336770,
                        ('Fe2O3', 'f'): -843.1471000,
                        ('Fe2O3', 'g'): 228.3548000,
                        ('Fe2O3', 'h'): -825.5032000,
                        ('CH4', 'a'): -0.7030290,
                        ('CH4', 'b'): 108.4773000,
                        ('CH4', 'c'): -42.5215700,
                        ('CH4', 'd'): 5.8627880,
                        ('CH4', 'e'): 0.6785650,
                        ('CH4', 'f'): -76.8437600,
                        ('CH4', 'g'): 158.7163000,
                        ('CH4', 'h'): -74.8731000,
                        ('CO2', 'a'): 24.9973500,
                        ('CO2', 'b'): 55.1869600,
                        ('CO2', 'c'): -33.6913700,
                        ('CO2', 'd'): 7.9483870,
                        ('CO2', 'e'): -0.1366380,
                        ('CO2', 'f'): -403.6075000,
                        ('CO2', 'g'): 228.2431000,
                        ('CO2', 'h'): -393.5224000,
                        ('H2O', 'a'): 30.0920000,
                        ('H2O', 'b'): 6.8325140,
                        ('H2O', 'c'): 6.7934350,
                        ('H2O', 'd'): -2.5344800,
                        ('H2O', 'e'): 0.0821390,
                        ('H2O', 'f'): -250.8810000,
                        ('H2O', 'g'): 223.3967000,
                        ('H2O', 'h'): -241.8264000
                                                }                        
    # solid phase only physical xteristics - most acquired from Abad et al. 2007
        self.radg = Param(default = 2.6e-7,
                          doc = 'rep. particle grain radius within OC part, m')
        self.dp = Param(default = 1.5e-3,
                        doc = 'Diameter of solid particles, m')

        self.rho_sol = Param(default = 3251.75, 
                             doc = 'Bulk density of solid particles, kg/m^3') 
        self.k_sol = Param(default = 12.3,
                             doc = 'Thermal conductivity of solid particles, J/mKs') #EPAT
        self.rhom = Param(default = 32811,
                          doc = 'molar density of carrier particle, mol/m^3')
        ####

        self.k0 = Param(self.rxn_idx, default = 8e-4,
                        doc = 'pre-exponential factor, mol^(1-N_rxn)m^(3*N_rxn -2)/s')
        # ^ is this a factor in the reaction rate ?
        # A: it is a factor in the rate constant expression, which itself
        # is a factor in the reaction rate 'rgen'
        self.E = Param(self.rxn_idx, default = 4.9e4,
                       doc = 'Activation energy, kJ/mol')
        # 
        self.N_rxn = Param(self.rxn_idx, default = 1.3,
                           doc = 'reaction order in gas species, (-)')
        # ^ how is this determined/why is it significant ?
        # A: experimentally, from Abad et al, I assume
        self.b_rxn = Param(self.rxn_idx, default = 12,
                           doc = 'reaction stoich coeff, (-)')     
        # Coefficient of Fe2O3, significant because this is the 'reference component'
        # for sake of calculating reaction rate
        # Available volume for reaction - from EPAT report (1-ep)'
        self.a_vol = Param(within = Reals, default = 0.28, 
                         doc= 'available reaction vol. per vol. of OC')
        # ^ must be determined experimentally... based on what theory?
        #Viscosity constants - Reference: Perry and Green Handbook; McGraw Hill, 2008
        self.mu_param = {('CH4','a'):5.2546e-7,('CH4','b'):0.59006,('CH4','c'):105.67,('CH4','d'):0,
                         ('CO2','a'):2.148e-6,('CO2','b'):0.46,('CO2','c'):290,('CO2','d'):0,
                         ('H2O','a'):1.7096e-8,('H2O','b'):1.1146,('H2O','c'):0,('H2O','d'):0}
        # What equation are these plugged into? valid over what range?
        
        #Thermal conductivity constants - Reference: Perry and Green Handbook; McGraw Hill, 2008
        self.k_param = {('CH4','a'):8.3983e-6,('CH4','b'):1.4268,('CH4','c'):-49.654,('CH4','d'):0,
                        ('CO2','a'):3.69,('CO2','b'):-0.3838,('CO2','c'):964,('CO2','d'):1.86e6,
                        ('H2O','a'):6.204e-6,('H2O','b'):1.3973,('H2O','c'):0,('H2O','d'):0}
        

    def _make_params(self):        
        """
        Create the paramters and sets of the MB model.
        """
        # Declare Immutable Parameters
        self.pi = Param(default=2*acos(0), doc='pi')
        self.R  = Param(default=8.314459848,
                      doc='Universal Gas Constant [J/mol/K]') 
        self.g = Param(default=9.80665,
                       doc='Gravitational Acceleration [m/s2]')       
        # Heat capacities and other parameters for the energy balances
        # These are computed using the inlet conditions, in a separate 
        # Matlab file
#        self.Cpg = Param(default=2.24301,
#                         doc='Gas heat capacity [kJ/kg/K]')

# parameters for transfer across boundaries ##
        self.hf_in = Param(default=0.42757,
                        doc='Gas-solid heat transfer coeff. [kJ/(m2Ks)]')
        self.hw_in = Param(default=0.0857,
                doc='Internal convective heat transfer coeff., [kJ/(m2Ks)]')
        # ^ assume constant...
        
        self.Cpw = Param(default=0.5,
                         doc='Heat capacity of steel wall [kJ/kg/K]')
        self.rhow = Param(default=8050,
                doc='Reactor wall denisty [kg/m3]')
        
        # h0, aw, aw1 and Tamb not used when adiabatic
        self.h0 = Param(default=0.28094736,
                doc='Wall-ambiance heat transfer coefficient, [kJ/(m2 K s)]')   
        self.aw = Param(default=3.2129, doc='[1/m]') # inner area of wall?
        self.aw1 = Param(default=3.2727, doc='[1/m]') # outer area of wall?
        self.Tamb = Param(default=298.15, doc='Ambient temperature [K]')
        # ^ in theory, could be time dependent (a disturbance vairble)
        # ... for another day ...
        
        self.w_th = Param(default=0.005,
                          doc='Reactor wall thickness [m]')
        self.k_steel = Param(
                default=16.3*1e-3, 
                doc='Thermal conductivity of the steel reactor wall [kW/m/K]')
        self.k_refractory = Param(
                default=1.2*1e-3, 
                doc='Thermal conductivity of the refractory lining [kW/m/K]')
        self.refractory_th = Param(
                default=0.2,
                doc='Refractory thickness [m]') 
###############################################

# These could all be time indexed mutable parameters, if air conditions
# are to be treated as disturbances...
        self.Cp_air = Param(
                default=1.005,
                doc='Heat capacity of air at the film temp. [kJ/kg/K]')
        self.mu_air = Param(
                default=1.96*1e-5,
                doc='Dynamic viscosity of air at the film temp. [Pa s]')
        self.k_air = Param(
                default=0.0271*1e-3,
                doc='Thermal conductivity of air at the film temp. [kW/m/K]')
        self.beta = Param(
                default=3.27*1e-3,
                doc='Thermal expansion coefficient of air at film temp. [1/K]')
        self.nu_air = Param(
                default=16.97*1e-6,
                doc='Kinetic viscosity of air at the film temp. [m2/s]') 
        #   /alfa/
        self.alfa_air = Param(
                default=22.873433*1e-6,
                doc='Thermal diffusivity of air at the film temp. [m2/s]')       
###
        self.emf = Param(default=0.45, mutable = True,
                         doc='Minimum fluidization porosity [-]') 
        
        self.tuning_param = Param(default=1, mutable=True)
        self.tuning_param2 = Param(default=1, mutable=True)

# these are in a separate function from property parameters because... ?
# they have nothing to do with the components of the model?
        
        
    def _make_domain(self):
        """
        Create the axial dimensionless domain for the model 
        (along length of reactor)
        """
        self.z = ContinuousSet(bounds=(0.0,1.0), domain=Reals,
                               initialize=self.fe_set,
                               doc='Axial dimensionless domain')      
        self.t = ContinuousSet(bounds=(0.0,self.horizon),domain=NonNegativeReals, doc='Time domain')
        # ^ should be supplied by top-level IDAES flowsheet block

        # will need a _make_time_domain function? Or is this handled by IDAES - I think so...

        
    def _make_vars(self):
        """
        Create the model variables. 
        """
        # Discretisation point location and size
        # independent of time
        self.l = Var(self.z, domain=Reals,
                     doc='Location of Discrete Elements [m]')  
        # ^ how are these variables calculated...?

        # Reactor bed dimensions and characteristics
        self.L = Var(
                domain=Reals, bounds=(0.0,100.0), initialize=1.5,
                doc='Reactor bed length [m]')
        self.Dr = Var(
                domain=Reals, bounds=(0.0,40.0), initialize=2,
                doc='Reactor inner diameter [m]')
        self.A_bed = Var(
                domain=Reals, bounds=(0.0,100.0), initialize=1,
                doc = 'Reactor cross-sectional area [m2]')
        # ^ does not get fixed, but gets specified by a constraint
        self.eps = Var(
                domain=Reals, bounds=(0.0,1.0), initialize=0.8,
                doc='Bed voidage [-]')                  
        # ^ why are these vars and not mutable params?
        # will be fixd for model solution, I imagine...
        
        # Gas Inlet Conditions   
        # In general, time varying
        # convention will be that t is the last index listed

        # these inlet conditions are the inputs to the process
        self.t.display()
        self.Gas_In_F = Var(self.t,
                domain=Reals, bounds=(0.0,1e5), initialize=100.0,
                    doc='Inlet gas molar flow rate [mol/s]')
        self.Gas_In_P = Var(self.t,
                domain=Reals, bounds=(0.0,1e2), initialize=1.05,
                    doc='Pressure of gas at gas inlet [bar]')
        self.Gas_In_Tg = Var(self.t,
                domain=Reals, bounds=(0.0,1000.0), initialize=298.15,
                    doc='Temperature of gas at gas inlet [K]')
        self.Gas_In_y = Var(
                self.GasList, 
                self.t,
                domain=Reals, bounds=(0.0,1.0), initialize=0.1,
                    doc='Mole fractions of gas species at gas inlet [mol/mol]')
        # again, these will (probably) not be allowed to vary during a solve
        # ^ fixed during setInputs function

        # Gas Outlet Conditions
        self.Gas_Out_P = Var(
                self.t,
                domain=Reals, bounds=(1.0,1e2), initialize=1.05,
                    doc='Pressure of gas at gas outlet [bar]')
        
        # Solid Inlet Conditions   
        
        self.Solid_In_M = Var(
                self.t,
                domain=Reals, bounds=(0.0,1e5), initialize=2280,
                    doc='Inlet solid mass flow rate [kg/s]')
        self.Solid_In_Ts = Var(
                self.t,
                domain=Reals, bounds=(0.0,2000.0), initialize=293.15,
                    doc='Temperature of gas at gas inlet [K]')
        self.Solid_In_x =  Var(
                self.SolidList,
                self.t,
                domain=Reals, bounds=(0.0,1.0), initialize=0.1,
                doc='Mass fraction of solid species at solid inlet [kg/kg]')  
        # ^ fixed during setInputs function
        
        # Solid Outlet Conditions   
        self.Solid_Out_M = Var(
                self.t,
                domain=Reals, bounds=(0.0,1e5), initialize=2280,
                    doc='Outlet solid mass flow rate [kg/s]')
        self.Solid_Out_M_Comp = Var(
                self.SolidList, 
                self.t,
                domain=Reals, bounds=(0.0,1e5), 
                initialize=2280,		
                doc='Outlet solid mass component flow rate [kg/s]')         
        self.Solid_Out_Ts = Var(
                self.t,
                domain=Reals, bounds=(0.0,2000.0), initialize=293.15,
                    doc='Temperature of gas at gas outlet [K]')
        self.Solid_Out_x =  Var(
                self.SolidList, 
                self.t,
                domain=Reals, bounds=(0.0,1.0), initialize=0.1,
                doc='Mass fraction of solid species at solid outlet [kg/kg]') 
        # ^ fixed some 'doc' strings here
        
        # Gas Phase Concentrations
        self.CgT = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=1.0, bounds=(0.0,1e5),
                doc='Bulk gas total concentration [mol/m3]')
        self.Cg = Var(
                self.z, 
                self.GasList, 
                self.t,
                domain=NonNegativeReals, initialize=1.0,
                bounds=(0.0,1e5), 
                doc='Bulk gas component concentrations [mol/m3]')
        self.Ctrans = Var(
                self.z, 
                self.GasList, 
                self.t,
                domain=Reals, initialize=0.0, 
                doc='Moles transferred from Gas to Solid [mol/m3/s]')
        # ^ why is this named Ctrans? would Mtrans be more intuitive?
        self.G_flux = Var(
                self.z, 
                self.GasList, 
                self.t,
                domain=Reals, initialize=1,
                doc="Superficial gas flux, \
                    superfiial gas (velocity)*(component conc.) (mol/m2/s)")  
        # flux in z-direction along the column, right?

        self.X_gas = Var(
                self.t,
                domain=Reals, bounds=(0.0,1.0), initialize=0.0,
                doc="Gaseous fuel conversion [-]")
        # ^ should this be indexed by z? 
        # equation document given to me suggests no

        # Gas Mole Fractions
        self.y = Var(
                self.z, 
                self.GasList, 
                self.t,
                domain=Reals, #bounds=(0.0,1.0),
                initialize=0.1, doc='Gas phase mole fractions [-]')
        self.ytot = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=0.1, #bounds=(0.0,1.0),
                doc='Gas phase total mole fraction [-]')  
        
        # Gas flowrates
        self.Ftotal = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=100.0, bounds=(0.0,1e5),
                doc='Total molar flow rate of gas [mol/s]')
        self.F = Var(
                self.z, 
                self.GasList, 
                self.t,
                domain=Reals, initialize=100.0,
                doc='Molar flow rate of gas components [mol/s]')      
        self.Gas_M = Var(
                self.z, 
                self.GasList, 
                self.t,
                domain=Reals, initialize=100.0,
                doc='Mass flow rate of gas components [kg/s]')
        
        
        # Solid Phase Loading     
        self.qT = Var(
                    self.z, 
                    self.t,
                    domain=Reals, initialize=1.0, bounds=(0.0,1e5),
                    doc='Total solid loading [kg/m3]')
        # ^ what is this exactly? different from density?
        self.q = Var(
                    self.z, 
                    self.SolidList, 
                    self.t,
                    domain=Reals, initialize=1.0,
                    bounds=(0.0,1e5), doc='Solid loading [kg/m3]')
        self.qtrans = Var(
                    self.z, 
                    self.SolidList, 
                    self.t,
                    domain=Reals,initialize=0.0,
                    doc='Moles transferred from Gas to Solid [kg/m3/s]')
        self.S_flux = Var(
                self.z, 
                self.SolidList, 
                self.t,
                domain=Reals, initialize=1,
                doc="Superficial solids flux, \
                    superficial solids (velocity)*(loading) (kg/m2/s)")  
        self.X_OC = Var(
                self.t,
                domain=Reals, bounds=(0.0,1.0), initialize=0.0,
                doc="Oxygen carrier conversion [-]")
        
        # Solid mass fractions
        self.x = Var(
                self.z, 
                self.SolidList, 
                self.t,
                domain=Reals, #bounds=(0.0,1.0),
                initialize=0.1, doc='Solid phase mass fractions [-]')
        self.xtot = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=0.1, #bounds=(0.0,1.0),
                doc='Solid phase total mass fraction [-]')            
                
        # Mass flowrates
        self.Solid_M_total = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='Total mass flow rate of solids [kg/s]')
        self.Solid_M = Var(
                self.z, 
                self.SolidList, 
                self.t,
                domain=Reals, 
                bounds=(0.0,1e7), initialize=1.0,
                doc='Mass flow rate of solid components [kg/s]') 
        self.Solid_F_total = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,1e7), initialize=1.0,
                doc='Total molar flow rate of solids [mol/s]')
        self.Solid_F = Var(
                self.z, 
                self.SolidList, 
                self.t,
                domain=Reals, 
                bounds=(0.0,1e7), initialize=1.0,
                doc='Molar flow rate of solid components [mol/s]') 
        self.mFe_mAl = Var(
                self.z,
                self.t,
                domain=Reals, initialize=1,
                doc='Fe to Al mass ratio [-]')
        
        # Temperatures and Pressures
        self.P = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,1e2), initialize=1.05,
                doc='Pressure [bar]')        
        self.Tg = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,2000.0), initialize=298.15,
                doc='Gas phase temperature [K]')
        self.Ts = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,2000.0), initialize=0.29815,
                doc='Solid phase temperature [kK]')           
        self.Tw = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,10000.0), initialize=298.15,
                doc='Wall temperature [K]')
        self.Tg_GS = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=0.0,
                doc = 'Gas to solid convective heat transfer [kJ/(m2.s)]')
        self.Tg_GW = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=0.0,
                doc = 'Gas to wall convective heat transfer [kJ/(m2.s)]')  
        self.Tg_refractory = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=0.0,
                doc = 'Heat transfer through refractory wall [kJ/(m2.s)]')
        self.Ts_dHr= Var(
                self.z, 
                self.t,
                domain=Reals, initialize=0.0,
                doc = 'Heat transfer due to heat of reaction [kJ/(m2.s)]')
        self.Tw_GW = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=0.0,
                doc = 'Heat transfer through reactor wall [kJ/(m2.s)]') 
        self.Tw_Wamb = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=0.0,
                doc = 'Heat transfer from reactor wall to ambient [kJ/(m2.s)]') 
        # ^ why are these heat transfer (heat flux) variables named 'T*'?
        
        # Velocities
        self.vg_in = Var(
                self.t,
                  domain=Reals, bounds=(-0.001,1e3), initialize=0.05,
                  doc='Inlet superficial gas velocity [m/s]')  # at inlet
        # ^ indexed by time; this is a disturbance/input to the process
        # (^ or is this an algebraic variable entirely determined by 
        # inlet conditions?)
        self.vg = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(-0.001,1e3), initialize=0.05,
                doc='Superficial gas velocity [m/s]')
        self.umf = Var(
                self.z, 
                self.t,
                domain=Reals,bounds=(-0.001,1e3), initialize=1.00584,
                doc='Minimum fluidization velocity [m/s]')
        self.vs = Var(
                self.t,
                  domain=Reals, initialize=5.00,
                  doc='Superficial solids velocity (downwards) [mm/s]') # not negative
        # ^ assumed constant along reactor... 
        # not constant with time; an input/disturbance
        self.v_diff=Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,1e5),
                doc='Variable to store "umf-vg", used in optimization')                
        
        # Dimensionless numbers
        self.Rep = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,1e8), initialize=1.0,
                doc='Particle Reynolds number (dimensionless)')
        self.Re = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,1e8), initialize=1.0,
                doc='Bed Reynolds number (dimensionless)')
        self.Pr = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='Prandtl number (dimensionless)')
        self.Pr_ext = Var(
                domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='External (air) Prandtl number (dimensionless)')
        # ^ depends only on ambient conditions, so independent of 
        # z,t
        self.Ra = Var(
                self.t,
                domain=Reals, bounds=(0.0,1e17), initialize=1.0,
                doc='Rayleigh number (dimensionless)')        
        # ^ depends on ambient conditions and wall temperature
        # at top of reactor, so dependent on t, but not z
        self.Nu = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='Nusselt number (dimensionless)')
        self.Nuw = Var(
                self.z, 
                self.t,
                domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='Nusselt number of wall (dimensionless)')
        self.Nu_ext = Var(
                self.t,
                domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='External Nusselt number (dimensionless)')
        self.hf = Var(
                self.z, 
                self.t,
                initialize=0.47502,
                        doc='Gas-solid heat transfer coeff. [kJ/(m2Ks)]')
        self.hw = Var(
                self.z, 
                self.t,
                initialize=0.21343,
                doc='Internal convective heat transfer coeff., [kJ/(m2Ks)]')
        self.hext = Var(
                self.t,
                bounds=(0.0,1e5), initialize=0.0018,
                doc='External natural convective HTC [kJ/(m2Ks)]')
        self.hext2 = Var(
                self.z, 
                self.t,
                bounds=(0.0,1e5), initialize=0.02,
                doc='External natural convective HTC [kJ/(m2Ks)]')
        # ^ why are there two?
        # and why is the second indexed by time
        self.U = Var(
                self.z, 
                self.t,
                bounds=(0.0,1e5), initialize=0.0017,
                doc='Overall heat transfer coeff. (including refractory), \
                [kJ/(m2Ks)]')
        self.Uw = Var(
                self.z, 
                self.t,
                bounds=(0.0,1e5), initialize=0.0017,
                doc='Overall heat transfer coeff., [kJ/(m2Ks)]')

        # Gas phase thermal flux
        self.Gh_flux = Var(
                self.z, 
                self.t,
                domain=Reals, initialize=100,
                doc='Gas phase enthalpy flux [kJ/(m2.s)]') 

        self.Sh_flux = Var(
                self.z,
                self.t,
                domain=Reals,
                initialize=100,
                doc='Solid phase enthalpy flux [kJ/m2.s]')
              
        # Dummy variables
#        self.dummy = Var(
#                self.z, domain=Reals, initialize = 1.0, 
#                doc='Dummy variable used to store dP/dz')

    def _make_init_values(self):
        self.Cg_0 = Param(self.z,self.GasList,mutable=True,default=1.0)
        self.q_0 = Param(self.z,self.SolidList,mutable=True,default=1.0)
        self.Tg_0 = Param(self.z,mutable=True,default=298.15)
        self.Ts_0 = Param(self.z,mutable=True,default=0.29815)
        # this should go somewhere else... :
        # why doesn't this work???
        print(self)
        #self.diff_vars_t = [self.Cg,self.q,self.Tg,self.Ts]
    
    def _make_vars_props(self):
        """
        Create properties variables. These would belong to a property package, 
        but have been moved here because it somehow improved convergence. Will 
        be moved back out while transitioning to the new framework.
        """ 

        # assuming these variables all depend on thermodynamic state (T,P,comp.),
        # they should all be indexed by time.

#        self.H_comp_s = Var(self.z, self.comp, domain = Reals, initialize = 1.0, 
#                        doc = 'Solid phase component Enthalpy at Ts(K) - J/(mol comp)')
# ^ replaced by expression
        self.DH_rxn_s = Var(self.z, self.t, domain = Reals, initialize = 1e-3, 
                        doc = 'Heat of rxn at system T, kJ/(mol reaction)')
        # ^ specified, except at z = 0, by eq_p2

        # specified by p4
        self.cp_sol = Var(self.z, self.t, domain = Reals, initialize = 1.0,
                        doc = 'Heat capacity of solid particles, J/kg.K')

        # specified by p5
        self.MW_vap = Var(self.z, self.t, domain=Reals, initialize = 0.03,
                          doc = 'mol wt. of gas - kg/mol')
        
        self.rho_vap = Var(self.z, self.t, domain=Reals, initialize = 1.0,
                        doc = 'gas mass density - kg/m3')
        self.mu_vap = Var(self.z, self.t, domain=Reals, initialize = 1e-2,
                        doc = 'dynamic viscosity of gas, - mPa s = g/m/s')
        self.cp_gas = Var(self.z, self.t, domain=Reals, initialize = 1.0,
                        doc='Gas phase heat capacity [kJ/kg/K]') 
        self.cp_vap = Var(self.z, self.t, domain = Reals, initialize = 1.0,
                        doc = 'heat capacity of gas - J/mol/K')
        self.cv_vap = Var(self.z, self.t, domain = Reals, initialize = 1.0,
                        doc = 'Constant vol. heat capacity of gas - J/mol/K')
        self.k_cpcv = Var(self.z, self.t, domain = Reals, initialize = 1.0,
                        doc = 'specific heat ratio (-)')        
        self.k_vap = Var(self.z, self.t, domain = Reals, initialize = 1e-1, 
                        doc = 'thermal conductivity of gas, J/mKs')
        self.X = Var(self.z, self.t, domain=Reals, initialize=0.0, 
                        doc = 'fraction of metal oxide converted')
        self.X_term = Var(self.z, self.t, domain = Reals, initialize = 1.0, 
                        doc = 'reformulation term for X to help eqn scaling') 
        self.k = Var(self.z, self.rxn_idx, self.t, domain=Reals, initialize=1.0, 
                        doc = 'kinetic rate constant')
        self.r_gen = Var(self.z, self.rxn_idx, self.t, domain=Reals, initialize = 1.0, 
                        doc = 'gen. rate expression, units= mol_rxn/kgOC.s')          
        self.rg = Var(self.z, self.GasList, self.t, domain=Reals, initialize = 0.0,
                        doc = 'gcomp. total rate expression, units= mol/m^3.s')        
        self.rs = Var(self.z, self.SolidList, self.t, domain=Reals, initialize = 0.0,
                      doc = 'scomp. rate expression, units= kg/m^3.s')
        # ^ total rates of accumlation of gas and solid phases resspectively?
    
    
    def _make_derivatives(self):
        """
        Make derivative variables
        """
        self.dG_fluxdz = DerivativeVar(self.G_flux, wrt=self.z)
        self.dS_fluxdz = DerivativeVar(self.S_flux, wrt=self.z)

        # need derivatives of C(z,j,t), q(z,j,t), Tg(z,t), Ts(z,t)
        # wrt time

        self.dCgdt = DerivativeVar(self.Cg,wrt=self.t,initialize=0.0)
        self.dqdt = DerivativeVar(self.q,wrt=self.t,initialize=0.0)
        
        self.dPdz = DerivativeVar(self.P, wrt=self.z) 
        
#        self.dTgdz = DerivativeVar(self.Tg, wrt=self.z)
        #self.dTsdz = DerivativeVar(self.Ts, wrt=self.z)        
        self.dldz = DerivativeVar(self.l, wrt=self.z)
        # ^ what is the reason to define the derivative of length (in the z direction) 
        # with respect to z? ... the world may never know...

        
        self.dGh_fluxdz = DerivativeVar(self.Gh_flux, wrt=self.z)
        self.dSh_fluxdz = DerivativeVar(self.Sh_flux, wrt=self.z)

        self.dTsdt = DerivativeVar(self.Ts, wrt=self.t,initialize=0.0)
        self.dTgdt = DerivativeVar(self.Tg, wrt=self.t,initialize=0.0)


    def _make_constraints_props(self):
        """ 
        Create property constraints. These would belong to a property package, 
        but have been moved here because it somehow improved convergence. Will 
        be moved back out while transitioning to the new framework.
        """
        # any rule defining a temperature (or Pressure/comp)-
        # dependent parameter must be indexed by time

        # why is the model called 'b' in these constraints? 

        # why are some equalities defined as expressions and some as constraints?


        def rule_eq_p1(b,z,i,t):
            return 1e3*(b.cp_param[i,'a']*(1e3*b.Ts[z,t]/1000)
                        + b.cp_param[i,'b']*((1e3*b.Ts[z,t]/1000)**2)/2
                        + b.cp_param[i,'c']*((1e3*b.Ts[z,t]/1000)**3)/3
                        + b.cp_param[i,'d']*((1e3*b.Ts[z,t]/1000)**4)/4
                        - b.cp_param[i,'e']/(1e3*b.Ts[z,t]/1000) 
                        + b.cp_param[i,'f']
                        - b.cp_param[i,'h'])
        self.H_comp_s = Expression(self.z, self.comp, self.t, rule=rule_eq_p1,
                                   doc = 'Component enthalpy [J/mol]')      
        

        def rule_eq_p2(b,z,t):
            # if H_comp, dH_rxn are defined at z=0, why not define this constraint there?
            if z == b.z.first():
                return Constraint.Skip
            else:
                return (1e3*b.DH_rxn_s[z,t]) == sum(b.stoic[i]*(b.H_comp_s[z,i,t] 
                                        + b.Hf_comp[i]) for i in b.comp)
        self.eq_p2 = Constraint(self.z, self.t, rule = rule_eq_p2,
                                doc = 'Heat of reaction')     
        # ^ why is constraint not specified at z=0
        
        # Solid properties
        def rule_eq_p3(b,z,i,t):
            return b.cp_param[i,'a'] + b.cp_param[i,'b']*(1e3*b.Ts[z,t]/1000) \
                        + b.cp_param[i,'c']*(1e3*b.Ts[z,t]/1000)**2 \
                        + b.cp_param[i,'d']*(1e3*b.Ts[z,t]/1000)**3 \
                        + b.cp_param[i,'e']/((1e3*b.Ts[z,t]/1000)**2)
        self.cp_comp_s = Expression(self.z, self.comp, self.t, rule=rule_eq_p3,
                                    doc = 'Solid comp. heat capacity [J/(mol.K)')
        
        def rule_eq_p4(b,z,t):
            return b.cp_sol[z,t] == sum(b.cp_comp_s[z,j,t]*(1000/b.MW[j])
                                      *b.x[z,j,t] for j in self.SolidList)
        self.eq_p4 = Constraint(self.z, self.t, rule=rule_eq_p4,
                                doc = 'Solid phase heat capacity')  
        # specifies cp_sol
        
        # Gas mixture molecular weight
        def rule_eq_p5(b,z,t):
            return b.MW_vap[z,t] == 1e-3*sum(b.y[z,i,t]*b.MW[i] for i in b.GasList)
        self.eq_p5 = Constraint(self.z, self.t, rule=rule_eq_p5,
                                doc = 'Gas mixture molecular weight')
        # ^ specifies MW_vap

        # Gas phase density
        def rule_eq_p6(b,z,t):
            return b.rho_vap[z,t] == b.MW_vap[z,t]*b.P[z,t]/(b.R*1e-5*b.Tg[z,t])
        self.eq_p6 = Constraint(self.z, self.t, rule=rule_eq_p6,
                                doc = 'Gas phase density')
        # Components viscosity 
        def rule_eq_p7(b,z,i,t):
            return  b.mu_param[i,'a']*(b.Tg[z,t]**b.mu_param[i,'b']) \
                        / ((1 + (b.mu_param[i,'c']/b.Tg[z,t])) \
                        + (b.mu_param[i,'d']/(b.Tg[z,t]**2))) 
        self.mu_comp = Expression(self.z,self.GasList, self.t, rule=rule_eq_p7,
                                  doc = 'Gas comp. dynamic viscosity [kg/ms]')
        
        # Viscosity of gas mix
        def rule_eq_p8(b,z,t):
            return 1e6*(1e-3*b.mu_vap[z,t]) == 1e6*sum(b.y[z,i,t]*b.mu_comp[z,i,t] \
                                    /(sum(b.y[z,j,t]*(b.MW[j]/b.MW[i])**0.5 \
                                    for j in b.GasList)) for i in b.GasList)
        self.eq_p8 = Constraint(self.z, self.t, rule=rule_eq_p8,
                                 doc = 'Gas mixture dynamic viscosity')

        # Heat capacity
        def rule_eq_p9(b,z,i,t):
            return b.cp_param[i,'a'] + b.cp_param[i,'b']*(b.Tg[z,t]/1000) \
                        + b.cp_param[i,'c']*(b.Tg[z,t]/1000)**2 \
                        + b.cp_param[i,'d']*(b.Tg[z,t]/1000)**3 \
                        + b.cp_param[i,'e']/((b.Tg[z,t]/1000)**2)
        self.cp_comp_g = Expression(self.z, self.GasList, self.t, 
                rule=rule_eq_p9, doc = 'Gas comp. heat capacity [J/(mol.K)')
            
        def rule_eq_p10(b,z,t):
            return b.cp_vap[z,t] == sum(b.cp_comp_g[z,j,t]*b.y[z,j,t] 
                    for j in b.GasList)
        self.eq_p10 = Constraint(self.z, self.t, rule=rule_eq_p10,
                                 doc = 'Gas mixture heat capacity')         
        
        def rule_eq_p11(b,z,t):
            return b.cp_gas[z,t]*b.MW_vap[z,t] == b.cp_vap[z,t]*1e-3
        self.eq_p11 = Constraint(self.z, self.t, rule=rule_eq_p11,
                                  doc = 'Gas mixture heat capacity [J/(kg.K)]')
        
        def rule_eq_p12(b,z,t):
            return b.cv_vap[z,t] == b.cp_vap[z,t] - b.R
        self.eq_p12 = Constraint(self.z, self.t, rule=rule_eq_p12,
                                  doc = 'Constant vol. heat capacity of gas')

        def rule_eq_p13(b,z,t):
            return b.k_cpcv[z,t] ==  b.cp_vap[z,t]/b.cv_vap[z,t]
        self.eq_p13 = Constraint(self.z, self.t, rule=rule_eq_p13,
                                  doc = 'Specific heat ratio of gas')        
                    
        # Thermal conductivity of gas  
        def rule_eq_p14(b,z,i,t):
            return b.k_param[i,'a']*(b.Tg[z,t]**b.k_param[i,'b']) \
                        / ((1 + (b.k_param[i,'c']/b.Tg[z,t])) \
                        + (b.k_param[i,'d']/(b.Tg[z,t]**2)))                         
        self.k_comp = Expression(self.z, self.GasList, self.t, rule=rule_eq_p14,
                                 doc = 'Gas comp. thermal conductivity [J/(m.K.s)]')

        def rule_eq_p15(b,z,i,j,t):
            return (1 + ((b.k_comp[z,j,t]/b.k_comp[z,i,t])**0.5) \
                    # ^ probably a typo..., should be different components
                        * ((b.MW[j]/b.MW[i])**0.25))**2 \
                        / (8*(1+(b.MW[j]/b.MW[i])))**0.5    
        self.A_bin = Expression(self.z, self.GasList, self.gcomp, self.t, 
                rule=rule_eq_p15, doc = 'Binary interaction param. btw gas species [-]')
            
        def rule_eq_p16(b,z,t):
            return self.k_vap[z,t] == sum(b.y[z,i,t]*b.k_comp[z,i,t] \
                                    /(sum(b.y[z,j,t]*b.A_bin[z,i,j,t]**0.5 \
                                    for j in b.GasList)) for i in b.GasList)
        self.eq_p16 = Constraint(self.z, self.t, rule=rule_eq_p16,
                                doc = 'Gas mixture thermal conductivity')    
        
        # rate expression constraints - overall rate and species specific rates
        def rule_eq_r1(b, z, i, t):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.k[z,i,t] == b.k0[i]*exp(-b.E[i]/(b.R*1e3*b.Ts[z,t]))
        self.eq_r1 = Constraint(self.z, self.rxn_idx, self.t, rule=rule_eq_r1, 
                                doc = 'kinetic rate constant eqn') 
        # ^ why is this equation skipped at z = 0, but k is defined at z = 0... ?
            
        """
        This equation is for calculating the conv. of fe2 to fe3 in the reactor. It is gotten
        by converting the EPAT conv. equation (based on mole flow) to a mass fraction equivalent.
        
            X =                         x_fe3
                ---------------------------------------------------------
                 x_fe3 + x_fe2*((stoich_fe3/-stoich_fe2)*(MW_fe3/MW_fe2))
                
        where  x_fe2,x_fe3 = wt. fraction of fe2 and fe3 at any
                    given point in the reactor (state variables)                     
        """  

        def rule_eq_r2(b,z,t):
            return 1e6*b.X[z,t]*(b.x[z,'Fe3O4',t] + (b.MW['Fe3O4']/b.MW['Fe2O3']) \
                        *(b.stoic['Fe3O4']/-b.stoic['Fe2O3'])*b.x[z,'Fe2O3',t]) \
                        == 1e6*b.x[z,'Fe3O4',t]
        self.eq_r2 = Constraint(self.z,self.t,rule=rule_eq_r2, 
                                doc = 'conversion of metal oxide eqn')  
        # ^ for some reason this is specified at z=0...

        def rule_eq_r3(b,z,t):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.X_term[z,t]**3 == (1-b.X[z,t])**2   
        self.eq_r3 = Constraint(self.z, self.t, rule = rule_eq_r3)        
        # X_term is defined for all z and never fixed, so why is this constraint skipped...
      
        
        # General rate expression for Fe2O3 reacting with CH4
        # why skip all reaction rates at z = 0
        # no reaction happens there, I assume, but all the variables are defined
        # something to do with the physics/geometry, or with the equations?
        def rule_eq_r4(b,z,i,t):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return (1e-3*b.r_gen[z,i,t])*1e6== 1e6*b.scale*b.x[z,'Fe2O3',t]*(b.a_vol/b.MW['Fe2O3'])\
                                    *3*b.b_rxn[i]*b.k[z,i,t]*(b.Cg[z,'CH4',t]**b.N_rxn[i])\
                                    *b.X_term[z,t]/(b.rhom*b.radg)/(-b.stoic['Fe2O3'])  
        self.eq_r4 = Constraint(self.z, self.rxn_idx, self.t, rule=rule_eq_r4, 
                        doc = 'general rate expression (mole extent basis), \
                                units are mol_reaction/gOC.s')       
        # ^ specifies r_gen

        def rule_eq_r5(b,z,i,t):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.rs[z,i,t] == b.rho_sol*b.MW[i]*b.stoic[i]*(1e-3*b.r_gen[z,1,t])
        self.eq_r5 = Constraint(self.z, self.SolidList, self.t, rule=rule_eq_r5, 
                                    doc = 'comp specific rate expression, \
                                    units are kgOC/m3/s')
        # ^ specifies rs

        def rule_eq_r6(b,z,i,t):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.rg[z,i,t] == -1e3*b.rho_sol*b.stoic[i]*(1e-3*b.r_gen[z,1,t])
        self.eq_r6 = Constraint(self.z, self.GasList, self.t, rule=rule_eq_r6, 
                                    doc = 'comp specific rate expression for a rxn,\
                                    units are mol/m3/s')
        # ^ specifies rg
            
    #==========================================================================
    def _make_bed_geometry(self): 
        """
        Reactor Sizing and Design
        """
        # should not be indexed by time 

        # a1 - Finite element location
        def rule_eq_a1(b, z):
            if z == b.z.first():
                return b.l[z] == 0
            else:
                return b.dldz[z] == b.L
        self.eq_a1 = Constraint(self.z, rule=rule_eq_a1,
                                doc = 'Finite element location')
        # ^ specifies l 

        # a2 - Area of bed
        def rule_eq_a2(b):
            return b.A_bed == (b.pi*(b.Dr/2)**2)
        self.eq_a2 = Constraint(rule=rule_eq_a2,
                                doc = 'Bed cross-sectional area')         
        # specifies A_bed
                
    #========================================================================== 
    def _mass_balance(self):
        """
        Add the mass balance constraints. The boundary conditions are defined 
        under '_make_bdry_conds'
        """          
        # Differential constraints


        # Bulk gas component mole balance
        def rule_eq_b1(b, z, j, t):    # j indexes GasList, z indexes dz
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.eps*b.L*b.dCgdt[z,j,t] \
                            == -b.dG_fluxdz[z,j,t] - \
                            (1-b.eps)*b.Ctrans[z,j,t]*b.L 
        self.eq_b1 = Constraint(self.z, self.GasList, self.t, rule=rule_eq_b1,
                                  doc = 'Bulk gas component mole balance')
        # ^ specifies dCgdt
        
        def rule_eq_b2(b, z, j, t): 
            if z == b.z.first():
                return Constraint.Skip
            else:
                return (1-b.eps)*b.L*b.dqdt[z,j,t] \
                        == b.dS_fluxdz[z,j,t] + (1-b.eps)*b.qtrans[z,j,t]*b.L 
        self.eq_b2 = Constraint(self.z, self.SolidList, self.t, rule=rule_eq_b2,
                                 doc = 'Solid component dynamic mass balance')
        # ^ specifies dqdt
        
        # Moles of gas reacting (transferred from Gas to Solid) [mol/m3/s]
        def rule_eq_b3(b, z, j, t):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 1e3*b.Ctrans[z,j,t] == 1e3*b.rg[z,j,t]
        self.eq_b3 = Constraint(self.z, self.GasList, self.t, 
                                    rule=rule_eq_b3,
                                    doc = 'Moles of gas reacted')      
        # ^ specifies Ctrans
        
        # Moles of solid reacting (transferred from Gas to Solid) [kg/m3/s]
        def rule_eq_b4(b, z, j, t):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 1e3*b.qtrans[z,j,t] == 1e3*b.rs[z,j,t]
        self.eq_b4 = Constraint(self.z, self.SolidList, self.t,
                                    rule=rule_eq_b4,
                                    doc = 'Moles of solid reacted')
        # ^ specifies qtrans
    
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Calculate the total and component GAS molar flow rates
        def rule_eq_c1(b,z,j,t):
            return b.F[z,j,t] == b.A_bed*b.G_flux[z,j,t]
        self.eq_c1 = Constraint(self.z, self.GasList, self.t, rule=rule_eq_c1,
                               doc = 'Gas component molar flow rate')
        # ^ specifies F
        
        def rule_eq_c2(b,z,j,t):
            return b.Gas_M[z,j,t] == b.F[z,j,t]*b.MW[j]*1e-3
        self.eq_c2 = Constraint(self.z, self.GasList, self.t, rule=rule_eq_c2,
                                   doc = 'Gas component mass flow rate')
        # ^ specifies Gas_M from F
        
        def rule_eq_c3(b,z,t):
            return b.Ftotal[z,t] == sum(b.F[z,j,t] for j in b.GasList)
        self.eq_c3 = Constraint(self.z, self.t, rule=rule_eq_c3,
                                    doc = 'Gas mixture total molar flow rate')
        # ^ specifies Ftotal from F

        # Calculate the gas components concentrations
        def rule_eq_c4(b,z,j,t):
            return b.Cg[z,j,t]*(b.A_bed*b.vg[z,t]) == b.F[z,j,t]  
        self.eq_c4 = Constraint(self.z, self.GasList, self.t, rule=rule_eq_c4,
                                doc = 'Gas component concentration')        
        # ^ specifies Cg
                
        # Bulk gas total concentration
        def rule_eq_c5(b, z, t):
            return b.CgT[z,t]*(b.A_bed*b.vg[z,t]) == b.Ftotal[z,t]
        self.eq_c5 = Constraint(self.z, self.t, rule=rule_eq_c5,
                                 doc = 'Gas mixture total concentration (in mol/m^3?)')
        # ^ specifies CgT
        
        # Calculate the gas component mole fractions 
        def rule_eq_c6(b, z, j, t):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.y[z,j,t]*b.CgT[z,t] == b.Cg[z,j,t]
        self.eq_c6 = Constraint(self.z, self.GasList, self.t, rule=rule_eq_c6,
                               doc = 'Gas component mole fraction')
        # ^ specifies y from compositions
        # y at z=0 will be specified directly from 
        
        def rule_eq_c7(b,z,t):
            return b.ytot[z,t] == sum(b.y[z,j,t] for j in b.GasList) 
        self.eq_c7 = Constraint(self.z, self.t, rule=rule_eq_c7,
                                  doc = 'Total mole fraction of gas') 
        # ^ specifies ytot from y
        
        def rule_eq_c8(b,t):
            return self.X_gas[t] == (1 - self.Gas_M[1,'CH4',t]/self.Gas_M[0,'CH4',t])
        #self.eq_c8 = Constraint(expr=self.X_gas == 
        #                            (1 - self.Gas_M[1,'CH4']/self.Gas_M[0,'CH4']),
        #                            doc = 'Conversion of gas fuel')
        self.eq_c8 = Constraint(self.t, rule=rule_eq_c8, 
                doc = 'Conversion of gas fuel')
        # ^ specifies X_gas from molar flow rates

        def rule_eq_c18(b,z,t):
            return b.CgT[z,t] == sum( b.Cg[z,j,t] for j in b.GasList )
        self.eq_c18 = Constraint(self.z, self.t, rule=rule_eq_c18)
        # if c18 is active, one of c3, c4, or c5 should be deactivated (probable c5)
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Calculate the total and component SOLID mass flow rates
        def rule_eq_c9(b,z,j,t):
            return b.Solid_M[z,j,t] == b.A_bed*b.S_flux[z,j,t]
        self.eq_c9 = Constraint(self.z,self.SolidList,self.t,rule=rule_eq_c9,
                                     doc = 'Solid component mass flow rate')
        # ^ specifies Solid_M from flux (differential)
        
        def rule_eq_c10(b,z,t):
            return b.Solid_M_total[z,t] == sum(b.Solid_M[z,j,t] \
                                              for j in b.SolidList)
        self.eq_c10 = Constraint(self.z, self.t, rule=rule_eq_c10,
                                           doc = 'Total mass flow rate of solid')
        # ^ specifies Solid_M_total 
        
        def rule_eq_c11(b,z,j,t):
            return b.Solid_F[z,j,t]*(b.MW[j]*1e-3) == b.Solid_M[z,j,t]
        self.eq_c11 = Constraint(self.z, self.SolidList, self.t, rule=rule_eq_c11,
                                     doc = 'Solid component mole flow rate')
        # specifies Solid_F from mass flow rate
        
        def rule_eq_c12(b,z,t):
            return b.Solid_F_total[z,t] == sum(b.Solid_F[z,j,t] \
                                              for j in b.SolidList)
        self.eq_c12 = Constraint(self.z, self.t, rule=rule_eq_c12,
                                           doc = 'Total mole flow rate of solid')
        # ^ specifies total molar flow rate
        
        # Calculate the solid components loading
        def rule_eq_c13(b,z,j,t):
            return b.q[z,j,t]*(1e-3*b.vs[t]) == b.S_flux[z,j,t]
        self.eq_c13 = Constraint(self.z, self.SolidList, self.t, rule=rule_eq_c13,
                               doc = 'Solid components loading')        
        # specifies q from vs & flux
        # or, alternatively, specifies S_flux from q and vs (constant)
        
        # Total solid loading [kg/m3]
        def rule_eq_c14(b, z, t):
            return b.qT[z,t] == sum(b.q[z,j,t] for j in b.SolidList) 
        self.eq_c14 = Constraint(self.z, self.t, rule=rule_eq_c14,
                                doc = 'Total solid loading')
        # ^ specifies total solid loading
        
        # Calculate the solid phase mass fractions 
        def rule_eq_c15(b, z, j, t):
            return b.x[z,j,t]*b.qT[z,t] == b.q[z,j,t]  
        self.eq_c15 = Constraint(self.z, self.SolidList, self.t, rule=rule_eq_c15,
                               doc = 'Solid component mass fractions')
        # ^ specifies x from densities
        
        def rule_eq_c16(b,z,t):
            return b.xtot[z,t] == sum(b.x[z,j,t] for j in b.SolidList) 
        self.eq_c16 = Constraint(self.z, self.t, rule=rule_eq_c16,
                                  doc = 'Total mass fraction of solid')  
        # ^ specifies xtot
        # is it not just 1? ... why is this not redundant...? 
        # exist components not in the solidlist?
        
        def rule_eq_c17(b,t):
            return b.X_OC[t] == 1 - self.Solid_M[0,'Fe2O3',t]/self.Solid_M[1,'Fe2O3',t]
        #self.eq_c17 = Constraint(expr=self.X_OC == 
        #                1 - self.Solid_M[0,'Fe2O3']/self.Solid_M[1,'Fe2O3'],
        #                doc = 'Oxygen carrier conversion')
        self.eq_c17 = Constraint(self.t, rule=rule_eq_c17,
                        doc = 'Oxygen carrier conversion')
        # ^ specifies X_OC 
                
    #==========================================================================     
    def _energy_balance(self):
        """ 
        Add the energy balance constraints. 
        """

        def rule_eq_d1(b, z, t):
            if z == b.z.first():
                return Constraint.Skip #The BC for Tg is under '_make_bdry_conds' 
            else:
                return b.eps*b.L*b.rho_vap[z,t]*b.cp_gas[z,t]*b.dTgdt[z,t] == \
                            - b.dGh_fluxdz[z,t] \
                            - b.Tg_GS[z,t]*b.L - b.Tg_GW[z,t]*b.L \
                            - b.Tg_refractory[z,t]*b.L
        self.eq_d1 = Constraint(self.z, self.t, rule=rule_eq_d1,
                                doc = 'Gas phase energy balance')
        # specifies dTgdt. 1e-6 is not a unit conversion...
        
        def rule_eq_d2(b, z, t):
            return b.Gh_flux[z,t] \
                   == b.rho_vap[z,t]*b.vg[z,t]*b.Tg[z,t]*b.cp_gas[z,t]
        self.eq_d2 = Constraint(self.z, self.t, rule=rule_eq_d2,
                                doc = 'Gas phase enthalpy flux') 
        # specifies Gh_flux from Tg (differential)
        
        def rule_eq_d3(b, z, t):
            return b.Tg_GS[z,t]*b.dp == b.tuning_param2 \
                              *6*(1-b.eps)*b.hf[z,t]*(b.Tg[z,t]-1e3*b.Ts[z,t])
        self.eq_d3 = Constraint(self.z, self.t, rule=rule_eq_d3,
                                   doc = 'Gas-solid convective heat transfer flux')
        # ^ specifies Tg_GS from gas/solid temps
        
        def rule_eq_d4(b, z, t):
            return b.Tg_GW[z,t]*b.Dr == b.hw[z,t]*(b.Tg[z,t] - b.Tw[z,t])*4
        self.eq_d4 = Constraint(self.z, self.t, rule=rule_eq_d4,
                                   doc = 'Gas-wall heat transfer flux')        
        # ^ specifies TgGW from gas/wall temps

        def rule_eq_d5(b, z, t):
            return b.Tg_refractory[z,t]*b.Dr == b.U[z,t]*(b.Tg[z,t] - b.Tw[z,t])*4
        self.eq_d5 = Constraint(self.z, self.t, rule=rule_eq_d5,
                                           doc = 'Heat flux through refractory')        
        # ^ specifies Tg_ref. from gas,wall temps...

        def rule_eq_d6(b, z, t):
            # shouldn't this be b.z.last()
            if z == b.z.first():
                return Constraint.Skip #The BC for Ts is under '_make_bdry_conds' 
            else:
                return (1-b.eps)*b.L*b.rho_sol*(b.cp_sol[z,t]*1e-3)*b.dTsdt[z,t] == \
                            b.dSh_fluxdz[z,t] \
                            + b.Tg_GS[z,t]*b.L \
                            + b.Ts_dHr[z,t]*b.L
        self.eq_d6 = Constraint(self.z, self.t, rule=rule_eq_d6,
                                doc = 'Solid phase energy balance')    
        # ... 1e-3 factor is to convert cp_sol from J to kJ /kgK 
               
        def rule_eq_d7(b, z, t):
            if z == b.z.first():
                return Constraint.Skip 
            else:
                return b.Ts_dHr[z,t] == b.tuning_param \
                                  *(1-b.eps)*b.rho_sol \
                                  *(-1e3*b.DH_rxn_s[z,t]) \
                                  *(1e-3*b.r_gen[z,1,t])
        self.eq_d7 = Constraint(self.z, self.t, rule=rule_eq_d7,
                                    doc = 'Heat of reaction flux')
        # ^ specifies TS_dHr 
        
        ########################
        def rule_eq_d8(b, z, t):
            return 0 == b.Tw_GW[z,t] - b.Tw_Wamb[z,t]
        self.eq_d8 = Constraint(self.z, self.t, rule=rule_eq_d8,
                                doc = 'Axial boundary condition for wall/ambient')
        def rule_eq_d9(b, z, t):
            return b.Tw_GW[z,t] == b.aw*b.hw[z,t]*(b.Tg[z,t]-b.Tw[z,t])
        self.eq_d9= Constraint(self.z, self.t, rule=rule_eq_d9,
                                   doc = 'Heat transfer from gas to wall')
        def rule_eq_d10(b, z, t):
            return b.Tw_Wamb[z,t] == b.aw1*b.Uw[z,t]*(b.Tw[z,t]-b.Tamb)
            #/(b.rhow*b.Cpw)
        self.eq_d10 = Constraint(self.z, self.t, rule=rule_eq_d10,
                                     doc = 'Heat transfer from wall to ambient')
        ########################
        # ^ these three constraints specify Tw_Wamb, Tw_GW, and Tw from Tg

        def rule_eq_d11(b,z,t):
            return b.Sh_flux[z,t] == 1e-3*b.qT[z,t]*b.cp_sol[z,t]*(1e-3*b.vs[t])*(1e3*b.Ts[z,t])
        self.eq_d11 = Constraint(self.z, self.t, rule=rule_eq_d11,
                                    doc = 'Solid phase enthalpy flux')
                
    #==========================================================================  
    def _pressure_balance(self):
        """ 
        Add the pressure balance constraints, Ergun equation.
        """
        
        # Ideal Gas Law
        def rule_eq_e1(b,z,t):
            if z == 0:
                return b.P[z,t] == b.Gas_In_P[t]  #+ b.ErgunRHS
            else:
                return b.P[z,t] == b.CgT[z,t]*(b.R*1e-5)*b.Tg[z,t]
        self.eq_e1 = Constraint(self.z, self.t, rule=rule_eq_e1,
                               doc = 'Pressure calculation from ideal gas law')
        # ^ specify P from inlet and/or ideal law


        # Pressure drop to superficial velocity correlation
        def rule_eq_e2(b, z, t):
            if b.press_drop == 'Ergun': # Ergun equation
                if z == b.z.first():
                    return b.vg[z,t] == b.vg_in[t]     # Inlet velocity
                else:
                    return -b.dPdz[z,t]*1e5/150 == (((1e-3*b.mu_vap[z,t])
                                    *(1-b.eps)**2*(b.vg[z,t] + (1e-3*b.vs[t])) 
                                    /(b.dp**2*b.eps**3)) 
                                    + 1/150*(1.75*b.rho_vap[z,t]*(1-b.eps) 
                                    *(b.vg[z,t] + (1e-3*b.vs[t]))*(b.vg[z,t] + (1e-3*b.vs[t])) 
                                    /(b.dp*b.eps**3)))*b.L

                    #################################################################
                    ## make the assumption here that (vg+vs) > 0 to allow symbolic ##
                    ## differentiation with simpy in the                           ##
                    ## calculate_variable_from_constraint() function.              ##
                    ## This should be regarded as temporary, as I may want to      ##
                    ## simulate later a case with negative velocities.             ##
                    ##           DO NOT FORGET ABOUT THIS!!!                       ##
                    #################################################################
                                                    ##
                                                    ##
                                                    ##
                                                    ##
                                                    ##
                                                    ##
                                                    ##
                                                    ##

                    #return -b.dPdz[z,t]*1e5 == ((150.0*b.mu_vap[z,t] 
                    #                *(1-b.eps)**2*(b.vg[z,t] + b.vs[t]) 
                    #                /(b.dp**2*b.eps**3)) 
                    #                + (1.75*b.rho_vap[z,t]*(1-b.eps) 
                    #                *abs(b.vg[z,t] + b.vs[t])*(b.vg[z,t] + b.vs[t]) 
                    #                /(b.dp*b.eps**3)))*b.L
            elif b.press_drop == 'SimplifiedP': # Mondino et al. (2017)
                if z == b.z.first():
                    return b.vg[z,t] == b.vg_in[t]     # Inlet velocity
                else:
                    return -b.dPdz[z,t]*1e5 == (
                                    (b.rho_sol - b.rho_vap[z,t])
                                    *0.2*b.vg[z,t]*b.L)
            else:
                raise Exception('press_drop method not recognised.')
        self.eq_e2 = Constraint(self.z, self.t, rule=rule_eq_e2,
                                    doc = 'Pressure drop to superfical gas \
                                        velocity correlation')
        # ^ specifies dPdz
        
        # Compute the minimum fluidized velocity 
        def rule_e3(b,z,t):
            return ((1e-3*b.mu_vap[z,t])**2)*((1.75/(b.emf)**3)*(b.dp*b.umf[z,t] 
                       *b.rho_vap[z,t]/(1e-3*b.mu_vap[z,t]))**2 \
                       + (150*(1-b.emf)/(b.emf)**3)*(b.dp \
                       *b.umf[z,t]*b.rho_vap[z,t]/(1e-3*b.mu_vap[z,t]))) == \
                       b.dp**3*b.rho_vap[z,t] \
                       *(b.rho_sol-b.rho_vap[z,t])*b.g
        self.eq_e3 = Constraint(self.z, self.t, rule=rule_e3,
                                         doc = 'minimum fluidization velocity')               
        # ^ specifies the minimum fluidized velocity

        # Store the difference between the minimum fluidization velocity and 
        # gas velocity, to use as either post-solve check or in optimization
        def rule_eq_e4(b,z,t):
            return b.v_diff[z,t] == b.umf[z,t] - b.vg[z,t]
        self.eq_e4 = Constraint(self.z, self.t, rule=rule_eq_e4,
                                        doc = 'Velocity post-check for moving \
                                            bed operating regime')
        # ^ specifies v_diff
                
    #==========================================================================
    def _make_bdry_conds(self):
        """
        Boundary conditions for balance equations.
        And inlet velocity of gas and solids.
        """
        # all should be indexed by time
        
        # BC for gas component mole balance
        def rule_eq_f1(b,j,t):
            return 1e2*b.G_flux[0,j,t]*b.A_bed \
                    == 1e2*b.Gas_In_F[t]*b.Gas_In_y[j,t]
        self.eq_f1 = Constraint(self.GasList, self.t, rule=rule_eq_f1,
                                     doc = 'Boundary condition for gas \
                                         component mole balance')
        # ^ specifies G_flux[0]

        def rule_eq_f11(b,j,t):
            return b.Gas_In_y[j,t] == b.y[0,j,t]
        self.eq_f11 = Constraint(self.GasList, self.t, rule=rule_eq_f11,
                                 doc = 'Direct specification of y at z=0')
        
        # BC for solid components mass balance
        def rule_eq_f2(b,j,t):
            return 1e2*b.S_flux[1,j,t]*b.A_bed \
                    == 1e2*b.Solid_In_M[t]*b.Solid_In_x[j,t] 
#             Vbed = (1-eps)*Abed*L, but L=1 here because of scaling
        self.eq_f2 = Constraint(self.SolidList, self.t, rule=rule_eq_f2,
                                    doc = 'Boundary condition for solid \
                                        component mass balance')  
        # ^ specifies S_flux[1]

        # BC for gas phase energy balance
        def rule_eq_f3(b,t):
            return b.Tg[0,t] == b.Gas_In_Tg[t] 
        #self.eq_f3 = Constraint(expr=self.Tg[0] == self.Gas_In_Tg,
        #                           doc = 'Boundary condition for gas phase \
        #                               energy balance')
        self.eq_f3 = Constraint(self.t, rule=rule_eq_f3,
                                   doc = 'Boundary condition for gas phase \
                                       energy balance')
        # ^ Specifies Tg[0]
        
        # BC for solid phase energy balance
        def rule_eq_f4(b,t):
            return (1e3*b.Ts[1,t]) == b.Solid_In_Ts[t]
        #self.eq_f4 = Constraint(expr=self.Ts[1] == self.Solid_In_Ts,
        #                           doc = 'Boundary condition for solid phase \
        #                               energy balance')
        self.eq_f4 = Constraint(self.t, rule=rule_eq_f4,
                                   doc = 'Boundary condition for solid phase \
                                       energy balance')
        # ^ specifies Ts[1]

        # ---------------------------------------------------------------------
        # Some inlet conditions

        # Inlet velocities    
        def rule_eq_f5(b,t):
            return b.vg_in[t]*b.A_bed*b.rho_vap[0,t] == b.Gas_In_F[t]*b.MW_vap[0,t]
        #self.eq_f5 = Constraint(expr=self.vg_in \
        #                *(self.A_bed*self.rho_vap[0]) \
        #                == self.Gas_In_F*self.MW_vap[0],
        #                doc = 'Inlet gas velocity')  
        self.eq_f5 = Constraint(self.t, rule=rule_eq_f5,
                        doc = 'Inlet gas velocity')  
        # ^ specifies vg_in

        def rule_eq_f6(b,t):
            return (1e-3*b.vs[t])*b.A_bed*b.rho_sol == b.Solid_In_M[t]
        #self.eq_f6 = Constraint(expr=self.vs*(self.A_bed \
        #                            *self.rho_sol) == self.Solid_In_M,
        #                        doc = 'Inlet solid velocity')      
        self.eq_f6 = Constraint(self.t, rule=rule_eq_f6,
                                doc = 'Inlet solid velocity')      
        # ^ specifies vs

        # ---------------------------------------------------------------------
        # Gas outlet conditions - used when connecting to other unit models
        def rule_eq_f7(b,t):
            return b.Gas_Out_P[t] == b.P[1,t]
        #self.eq_f7 = Constraint(expr=self.Gas_Out_P == self.P[1],
        #                               doc = 'Gas outlet pressure')
        self.eq_f7 = Constraint(self.t, rule=rule_eq_f7,
                                       doc = 'Gas outlet pressure')
        # ^ should have more gas inlet/outlet conditions, such as those for 
        # composition, see constraint for solids below...
        # where are these?
        # the variables they would specify simply don't exist...
        
        # Solids outlet conditions - used when connecting to other unit models
        def rule_eq_f8(b,t):
            return b.Solid_Out_M[t] == b.Solid_M_total[0,t]
        #self.eq_f8 = Constraint(expr= self.Solid_Out_M 
        #                                 == self.Solid_M_total[0],
        #                                 doc = 'Solid outlet total mass flow')    
        self.eq_f8 = Constraint(self.t, rule=rule_eq_f8, 
                                         doc = 'Solid outlet total mass flow')    
        def rule_eq_f9(b,t): 
            return b.Solid_Out_Ts[t] == (1e3*b.Ts[0,t])
        #self.eq_f9 = Constraint(expr= self.Solid_Out_Ts 
        #                                 == self.Ts[0],
        #                                 doc = 'Solid outlet temperature') 
        self.eq_f9 = Constraint(self.t, rule=rule_eq_f9, 
                                         doc = 'Solid outlet temperature') 
        def rule_eq_f10(b,j,t):
            return b.Solid_Out_x[j,t] == b.x[0,j,t]
        self.eq_f10 = Constraint(self.SolidList, self.t, rule=rule_eq_f10,
                                         doc = 'Solid outlet mass composition')
        
    #==========================================================================        

    def _make_init_conditions(self):
        # some skipped constraints because initial conditions 
        # should not be enforced at the inlets

        def rule_eq_h1(b,z,i):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Cg[z,i,0] == b.Cg_0[z,i]
        self.eq_h1 = Constraint(self.z,self.GasList,rule=rule_eq_h1)

        def rule_eq_h2(b,z,j):
            #if z == b.z.first():
            if z == b.z.last():
                return Constraint.Skip
            else:
                return b.q[z,j,0] == b.q_0[z,j]
        self.eq_h2 = Constraint(self.z,self.SolidList,rule=rule_eq_h2)

        def rule_eq_h3(b,z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Tg[z,0] == b.Tg_0[z]
        self.eq_h3 = Constraint(self.z,rule=rule_eq_h3)

        def rule_eq_h4(b,z):
            #if z == b.z.first():
            if z == b.z.last():
                return Constraint.Skip
            else:
                return (b.Ts[z,0]) == b.Ts_0[z]
        self.eq_h4 = Constraint(self.z,rule=rule_eq_h4)






    

    #==========================================================================        

    def _transfer_coefficients(self):
        """
        Calculate dimensionless numbers and heat transfer coefficients.
        """        
        # Particle Reynolds number   #*
        def rule_eq_g1(b, z, t):
            return b.Rep[z,t]*(1e-3*b.mu_vap[z,t]) == b.vg[z,t] \
                                *b.dp*b.rho_vap[z,t]
        self.eq_g1 = Constraint(self.z,self.t, rule=rule_eq_g1,
                                 doc = 'Particle Reynolds number')
        
        # Reynolds number
        def rule_eq_g2(b, z, t):
            return b.Re[z,t]*(1e-3*b.mu_vap[z,t]) == b.vg[z,t] \
                                *b.Dr*b.rho_vap[z,t]
        self.eq_g2 = Constraint(self.z, self.t, rule=rule_eq_g2,
                                doc = 'Bed Reynolds number')
        
        # Prandtl number - cp_gas [J/kg/K]
        def rule_eq_g3(b, z, t):
            return b.Pr[z,t] == (b.cp_gas[z,t]*1e3*
                                       (1e-3*b.mu_vap[z,t])/
                                       b.k_vap[z,t])
        self.eq_g3 = Constraint(self.z, self.t, rule=rule_eq_g3,
                                doc = 'Prandtl number of gas mixture in bed')

        # External Prandtl number (for air)
        self.eq_g4 = Constraint(expr=self.Pr_ext*self.k_air == 
                                    self.Cp_air*self.mu_air,
                                    doc = 'Prandtl number of ambient air')
        # ^ not indexed by time, but could be if ambient temperature 
        # is taken as a disturbance
        
        # External Rayleigh number (for air)
        # needs to be indexed by time because of dependence on 
        # wall temperature... other variables could be indexed by time
        # as disturbances...
        def rule_eq_g5(b,t):
            return b.Ra[t]*b.nu_air*b.alfa_air == b.g*b.beta*(b.Tw[1,t]-b.Tamb)*b.L**3
        #self.eq_g5 = Constraint(expr=self.Ra*self.nu_air*self.alfa_air == 
        #                self.g*self.beta*(self.Tw[1]-self.Tamb)*self.L**3.0,
        #                doc = 'Rayleigh number of ambient air')
        self.eq_g5 = Constraint(self.t, rule=rule_eq_g5, 
                        doc = 'Rayleigh number of ambient air')
        
         # Nusselt number of bed
        def rule_eq_g6(b, z, t):
            return b.Nu[z,t] == (2.0 + 1.1*(abs(b.Rep[z,t])**0.6)
                                   *(abs(b.Pr[z,t])**(1/3))) 
        self.eq_g6 = Constraint(self.z, self.t, rule=rule_eq_g6,
                                doc = 'Nusselt number of bed')
        
        # Nusselt number of wall
        def rule_eq_g7(b, z, t):
            return b.Nuw[z,t] == (12.5 + 0.048*b.Re[z,t])
        self.eq_g7 = Constraint(self.z, self.t, rule=rule_eq_g7,
                                 doc = 'Nusselt number of wall')

        # Exterior natural convection Nusselt number
        def rule_eq_g8(b,t):
            return b.Nu_ext[t] == 0.68 + (0.67*b.Ra[t]**(1/4)) /\
                                    ((1+(0.492/abs(self.Pr_ext))**(9/12))**(4/9))
        #self.eq_g8 = Constraint(
        #        expr=self.Nu_ext == 0.68 + (0.67*self.Ra**(1/4))
        #        /((1+(0.492/abs(self.Pr_ext))**(9/12))**(4/9)),
        #        doc = 'Exterior natural convection Nusselt number')
        self.eq_g8 = Constraint(self.t, rule=rule_eq_g8,
                doc = 'Exterior natural convection Nusselt number')
        
        # Gas-solid heat transfer coefficient
        def rule_eq_g9(b, z, t):
            return 1e6*b.hf[z,t]*b.dp == \
                                1e6*b.Nu[z,t]*b.k_vap[z,t]*1e-3
        self.eq_g9 = Constraint(self.z, self.t, rule=rule_eq_g9,
                                doc = 'Gas-solid heat transfer coefficient')
        
        # Gas-wall heat transfer coefficient
        def rule_eq_g10(b, z, t):
            return 1e6*b.hw[z,t]*b.eps*b.Dr == \
                            1e6*b.Nuw[z,t]*b.k_vap[z,t]*1e-3
        self.eq_g10 = Constraint(self.z, self.t, rule=rule_eq_g10,
                                doc = 'Gas-wall heat transfer coefficient')

        # Exterior natural convection Nusselt number
        def rule_eq_g11(b,t):
            return b.hext[t]*b.L == b.k_air*b.Nu_ext[t]
        #self.eq_g11 = Constraint(expr=self.hext*self.L == 
        #                          self.k_air*self.Nu_ext,
        #                          doc = 'External natural convective HTC')
        self.eq_g11 = Constraint(self.t, rule=rule_eq_g11, 
                                  doc = 'External natural convective HTC')
        # ^ must be indexed by time because Nu_ext is indexed by time 
        
        # Exterior natural convection Nusselt number 2
        def rule_eq_g12(b,z,t):
            return b.hext2[z,t]*b.L/b.nfe == b.k_air*b.Nu_ext[t]    
        self.eq_g12 = Constraint(self.z, self.t, rule=rule_eq_g12,
                                   doc = 'External natural convective HTC 2')
        # ^ why is this equation indexed by z? 
        # Will hext2 not be identical for all z?
        
        # Overall heat transfer coefficient (to use w/ refractory lining)
        def rule_eq_g13(b, z, t):
            return 1e8*b.U[z,t] == 1e8/((1/b.hw[z,t]) \
                        + b.refractory_th/b.k_refractory \
                        + b.w_th/b.k_steel \
                        + (1/b.hext2[z,t]))
        self.eq_g13 = Constraint(self.z, self.t, rule=rule_eq_g13,
                               doc = 'Overall HTC for refractory lining')

        # Overall heat transfer coefficient (wall-ambient)
        def rule_eq_g14(b, z, t):
            return 1e8*b.Uw[z,t] == 1e8/((1/b.hw[z,t]) \
                        + b.w_th/b.k_steel \
                        + (1/b.hext2[z,t]))
        self.eq_g14 = Constraint(self.z, self.t, rule=rule_eq_g14,
                                doc = 'Overall HTC for wall-to-ambient')

    #==========================================================================                
    def _dae_transform(self):
        """
        Apply DAE transformation to domain.
        """
        print('applying z-discretization')
        if self.dae_method == "OCLR":
            discretizer = TransformationFactory('dae.collocation')
            discretizer.apply_to(self, nfe=self.nfe, ncp=self.ncp,
                                 wrt=self.z, scheme='LAGRANGE-RADAU')          
        elif self.dae_method == "BFD1":
            discretizer = TransformationFactory('dae.finite_difference')
            discretizer.apply_to(self, nfe=self.nfe, wrt=self.z,
                                 scheme='BACKWARD')
        elif self.dae_method == "OCLL":
            discretizer = TransformationFactory('dae.collocation')
            discretizer.apply_to(self, nfe=self.nfe, ncp=self.ncp,
                                 wrt=self.z, scheme='LAGRANGE-LEGENDRE')
        else:
            raise Exception('DAE method type not recognised.')         
        print('nfe (z):\t',self.nfe)
        print('ncp (z):\t',self.ncp)

        '''
        Apply DAE transformation to time domain
        '''
        # should this be in a separate function? ...probably...
        # ... needs to be tested right now? ...

    def _t_dae_transform(self):
        t_discretizer = TransformationFactory('dae.collocation')
        print('applying t-discretization')
        print('nfe (t):\t',self.nfe_t)
        print('ncp (t):\t',self.ncp_t)
        t_discretizer.apply_to(self,
                nfe = self.nfe_t,
                ncp = self.ncp_t, 
                wrt = self.t,
                scheme = 'LAGRANGE-RADAU' )
        print('applied discretization')

        
            
    #==========================================================================            
'''
    def _initialize(blk, y=None, outlvl=0, optarg=None):
        # Initialisation routine for unit (default solver ipopt)
        
        # Create a solver and set solver options
        opt = SolverFactory('ipopt')
        
        # Strip bounds
        blk.strip_bounds()     
            
        if outlvl > 0:
            print("\n")
            print("----------------------------------------------------------")
            print("MB Fuel Reactor Initialization ...")
            print("\n")

        # Set Initial Values of State and Property Variables
        # Fix the temperatures, pressures and component fractions
        for z in blk.z:
            blk.P[z].fix(value(blk.Gas_In_P))
            blk.Tg[z].fix(value(blk.Solid_In_Ts)) #Because Thermal mass solid >> gas
            blk.Ts[z].fix(value(blk.Solid_In_Ts))    
            blk.Tw[z].fix(value(blk.Solid_In_Ts)) 
            for j in blk.GasList:
                blk.y[z,j].fix(value(blk.Gas_In_y[j]))
            for j in blk.SolidList:
                blk.x[z,j].fix(value(blk.Solid_In_x[j]))

        blk.Solid_Out_M.fix(value(blk.Solid_In_M))
        blk.Gas_Out_P.fix(value(blk.Gas_In_P))
        blk.Solid_Out_Ts.fix(value(blk.Solid_In_Ts))
        for j in blk.SolidList:
            blk.Solid_Out_x[j].fix(value(blk.Solid_In_x[j]))

        # ---------------------------------------------------------------------
        # Fix non-state derivative variables to ensure square problem
        blk.G_flux.fix(1.0)
        blk.S_flux.fix(1.0)
        blk.Gh_flux.fix(1.0)
                
        # ---------------------------------------------------------------------
        # Deactivate constraints
        # All property blocks (p and r constraints) are active        
        # All a constraints active
        
        blk.eq_b1.deactivate()
        blk.eq_b2.deactivate()
        blk.eq_b3.deactivate()
        blk.eq_b4.deactivate()
        
        blk.eq_c1.deactivate()
        blk.eq_c2.deactivate()
        blk.eq_c3.deactivate()
        blk.eq_c4.deactivate()
        blk.eq_c5.deactivate()
        blk.eq_c6.deactivate()
        blk.eq_c7.deactivate()
        blk.eq_c8.deactivate()
        blk.eq_c9.deactivate()
        blk.eq_c10.deactivate()
        blk.eq_c11.deactivate()
        blk.eq_c12.deactivate()
        blk.eq_c13.deactivate()
        blk.eq_c14.deactivate()
        blk.eq_c15.deactivate()
        blk.eq_c16.deactivate()
        blk.eq_c17.deactivate()
        
        blk.eq_d1.deactivate()
        blk.eq_d2.deactivate()
        blk.eq_d3.deactivate()
        blk.eq_d4.deactivate()
        blk.eq_d5.deactivate()
        blk.eq_d6.deactivate()
        blk.eq_d7.deactivate()
        blk.eq_d8.deactivate()
        blk.eq_d9.deactivate()
        blk.eq_d10.deactivate()
       
        blk.eq_e1.deactivate()
        blk.eq_e2.deactivate()
        blk.eq_e3.deactivate()
        blk.eq_e4.deactivate()
        
        blk.eq_f1.deactivate()
        blk.eq_f2.deactivate()
        blk.eq_f3.deactivate()
        blk.eq_f4.deactivate()
        blk.eq_f5.deactivate()
        blk.eq_f6.deactivate()
        blk.eq_f7.deactivate()
        blk.eq_f8.deactivate()
        blk.eq_f9.deactivate()
        blk.eq_f10.deactivate()
        
        blk.eq_g1.deactivate()
        blk.eq_g2.deactivate()
        blk.eq_g3.deactivate()
        blk.eq_g4.deactivate()
        blk.eq_g5.deactivate()
        blk.eq_g6.deactivate()
        blk.eq_g7.deactivate()
        blk.eq_g8.deactivate()
        blk.eq_g9.deactivate()
        blk.eq_g10.deactivate()
        blk.eq_g11.deactivate()
        blk.eq_g12.deactivate()
        blk.eq_g13.deactivate()
        blk.eq_g14.deactivate()

#%%        # ---------------------------------------------------------------------
        # 1st Initialisation Step - Properties Initialisation
        # Initialise solids mass fractions at equilibrium
        print()
        print("Initialise properties")
        for z in blk.z:                
#            if z != 0:
            blk.r_gen[z,1].fix(0.0)   
            blk.eq_r4.deactivate()

        # Solve mass balance equations with transfer off 
        ts = time.time()
        if outlvl > 0:
            print(blk,"- Step 1, initialize properties Time: ", end="")       

        results = opt.solve(blk,tee=True,symbolic_solver_labels=False,
                            keepfiles=False,options=optarg)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)    
            
#%%        # ---------------------------------------------------------------------
        # 2nd Initialisation Step - Hydrodynamics, mass & heat tranfer coeff
        # Initialise variables
        print()
        print("Initialise hydrodynamics") 
        for z in blk.z:
            blk.vg[z].fix(value(blk.Gas_In_F)*value(blk.MW_vap[0])
                          /(value(blk.eps)*(value(blk.pi)*(value(blk.Dr)/2)**2)
                          *value(blk.rho_vap[0])))
            blk.CgT[z].fix(value(blk.rho_vap[0])/value(blk.MW_vap[0])) 

        # Activate equations
        blk.eq_e1.activate()
        blk.eq_e3.activate()
        
        blk.eq_c5.activate()
        
        blk.eq_f5.activate()
        blk.eq_f6.activate()
        
        # unfix variables
        blk.P.unfix() 
        blk.umf.unfix()
        
        # Solve hydrodynamic, heat & mass transfer equations 
        ts = time.time()
        if outlvl > 0:
            print(blk,"- Step 2, initialize hydrodynamics Time: ", end="")       

        results = opt.solve(blk,tee=True,symbolic_solver_labels=False,
                            keepfiles=False,options=optarg)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)   


#%%        # ---------------------------------------------------------------------
        # 3rd Initialisation Step - Mass Balances
        # Initialise variables
        print()
        print("Initialize mass balances") 
                                    
        for z in blk.z:
            for j in blk.GasList:
                blk.G_flux[z,j] = value(blk.Gas_In_F)*value(blk.Gas_In_y[j]) \
                                    /(value(blk.pi)*(value(blk.Dr)/2)**2)
                blk.Ctrans[z,j] = value(blk.rg[z,j])
            for j in blk.SolidList:
                blk.S_flux[z,j] = value(blk.Solid_In_M)*value(blk.Solid_In_x[j]) \
                               /(value(blk.pi)*((value(blk.Dr)/2.0)**2)*1) 
                blk.qtrans[z,j] = value(blk.rs[z,j])
                #^ why the '==' here?
                # does this do anything?

        # Activate constraints - check the boundary conditions
        blk.eq_b1.activate()
        blk.eq_b2.activate()
        blk.eq_b3.activate()
        blk.eq_b4.activate()
        
        blk.eq_c1.activate()
        blk.eq_c2.activate()
        blk.eq_c3.activate()
        blk.eq_c4.activate()
        blk.eq_c6.activate()
        blk.eq_c7.activate()
        blk.eq_c8.activate()
        blk.eq_c9.activate()
        blk.eq_c10.activate()
        blk.eq_c11.activate()
        blk.eq_c12.activate()
        blk.eq_c13.activate()
        blk.eq_c14.activate()
        blk.eq_c15.activate()
        blk.eq_c16.activate()
        blk.eq_c17.activate()  
        
        blk.eq_e2.activate()
        blk.eq_e4.activate()

        blk.eq_f1.activate()
        blk.eq_f2.activate()
      
        # Unfix variables
        blk.vg.unfix()
        blk.CgT.unfix() 
        blk.G_flux.unfix()
        blk.S_flux.unfix()
        blk.y.unfix()
        blk.x.unfix()

        # Solve mass balance w/o rxn equations 
        ts = time.time()
        if outlvl > 0:
            print(blk,"- Step 3, initialize mass balances w/o rxn Time: ", end="")       

        results = opt.solve(blk,tee=True,symbolic_solver_labels=False,
                            keepfiles=False,options=optarg)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """) 

        # Activate reactions
        blk.eq_r4.activate() 
        
        # Unfix variables
        blk.r_gen.unfix()
        
        # Solve mass balance with rxn equations 
        ts = time.time()
        if outlvl > 0:
            print(blk,"- Step 3b, initialize mass balances with rxn Time: ", end="")       

        results = opt.solve(blk,tee=True,symbolic_solver_labels=False,
                            keepfiles=False,options=optarg)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)


#%%        # ---------------------------------------------------------------------
        # 4th Initialisation Step - Energy Balances and Transfer Eqns
        # Initialise variables
        print()
        print("Initialise energy balances") 
        
        # Initialize transfer variables
        for z in blk.z:
            blk.Rep[z] = value(blk.vg[z])*value(blk.dp) \
                        *value(blk.rho_vap[z])/value(blk.mu_vap[z])
            blk.Re[z] = value(blk.vg[z])*value(blk.Dr)*value(blk.rho_vap[z]) \
                        /value(blk.mu_vap[z])   
            blk.Pr[z] = 1e3*value(blk.cp_gas[z])*value(blk.mu_vap[z]) \
                        /value(blk.k_vap[z])
            blk.Nu[z] = 2.0 + 1.1*(abs(value(blk.Rep[z]))**0.6) \
                                   *(abs(value(blk.Pr[z]))**(1/3))
            blk.Nuw[z] = (12.5 + 0.048*value(blk.Re[z]))
            blk.hf[z] = 1e-3*value(blk.Nu[z])*value(blk.k_vap[z]) \
                            /value(blk.dp)
            blk.hw[z] = 1e-3*value(blk.Nuw[z])*value(blk.k_vap[z]) \
                        /(value(blk.eps)*value(blk.Dr))
            blk.hext2[z] = value(blk.nfe)*value(blk.k_air) \
                            *value(blk.Nu_ext)/value(blk.L)
            blk.U[z] = 1/((1/value(blk.hw[z])) \
                        + value(blk.refractory_th)/value(blk.k_refractory) \
                        + value(blk.w_th)/value(blk.k_steel) \
                        + (1/value(blk.hext2[z])))
            blk.Uw[z] = 1/((1/value(blk.hw[z])) \
                        + value(blk.w_th)/value(blk.k_steel) \
                        + (1/value(blk.hext2[z])))
            blk.Gh_flux[z] = value(blk.Ftotal[z])*value(blk.Tg[z]) \
                        *value(blk.cp_gas[z]) \
                        / ((value(blk.pi)*(value(blk.Dr)/2)**2))
        
        blk.Pr_ext = value(blk.Cp_air)*value(blk.mu_air)/value(blk.k_air)
        blk.Ra = value(blk.g)*value(blk.beta) \
                *(value(blk.Tw[1])-value(blk.Tamb))*value(blk.L)**3.0 \
                /(value(blk.nu_air)*value(blk.alfa_air))
        blk.Nu_ext = 0.68 + (0.67*value(blk.Ra)**(1/4)) \
                    /((1+(0.492/abs(value(blk.Pr_ext)))**(9/12))**(4/9))
        blk.hext = value(blk.k_air)*value(blk.Nu_ext)/value(blk.L)            

        # Activate energy balance eqns and heat transfer constraints
        blk.eq_d1.activate()
        blk.eq_d2.activate()
        blk.eq_d3.activate()
        blk.eq_d6.activate()
        blk.eq_d7.activate()
        
        blk.eq_f3.activate()
        blk.eq_f4.activate() 
        
        # Activate transfer constraints
        blk.eq_g1.activate()
        blk.eq_g3.activate()
        blk.eq_g6.activate()
        blk.eq_g9.activate()
        
        # Unfix variables
        blk.Tg.unfix()
        blk.Ts.unfix()    
        blk.Tw.unfix() 
        blk.Gh_flux.unfix()

        # fix wall related terms in gas energy balance to turn them off
        blk.Tg_refractory.fix(0.0)
        blk.Tg_GW.fix(0.0)
        
        # Solve energy balance equations w/o wall terms 
        ts = time.time()  
        if outlvl > 0:
            print(blk,"- Step 4, initialize energy bal w/o wall terms Time: ", end="")       

        results = opt.solve(blk,tee=True,symbolic_solver_labels=False,
                            keepfiles=False,options=optarg)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)        

#        # 4th Initialisation Step - Energy Balances and Transfer Eqns with wall terms
#        # Initialise variables
#        print()
#        print("Initialise energy balances and transfer with wall terms eqns")   
#
#        blk.eq_d4.activate()
#        blk.eq_d8.activate()
#        blk.eq_d9.activate()
#        blk.eq_d10.activate()
#            
#        blk.eq_g2.activate()
#        blk.eq_g4.activate()
#        blk.eq_g5.activate()
#        blk.eq_g7.activate()
#        blk.eq_g8.activate()
#        blk.eq_g10.activate()
#        blk.eq_g12.activate()
#        blk.eq_g14.activate()    
#
#        # Unfix wall terms
#        blk.Tg_GW.unfix()
#        
#        # Solve energy balance equations with wall terms 
#        ts = time.time()
#        if outlvl > 0:
#            print(blk,"- Step 4b, initialize energy bal with wall terms Time: ", end="")       
#
#        results = opt.solve(blk,tee=True,symbolic_solver_labels=False,
#                            keepfiles=False,options=optarg)
#        
#        if outlvl > 0:
#            print(value(time.time() - ts), " s")
#            print(results.solver.message)  
#            print("""    """)   
#
#
#        # 4th Initialisation Step - Energy Balances and Transfer Eqns with refractory terms
#        # Initialise variables
#        print()
#        print("Initialise energy balances and transfer with refractory terms eqns") 
#
#        # Activate refractory constraints
#        blk.eq_d5.activate()    
            
#        blk.eq_g13.activate()
#
#        # Unfix refractory terms
#        blk.Tg_refractory.unfix()
#        
#        # Solve energy balance equations with wall terms 
#        ts = time.time()
#        if outlvl > 0:
#            print(blk,"- Step 4c, initialize energy bal with refrac terms Time: ", end="")       
#
#        results = opt.solve(blk,tee=True,symbolic_solver_labels=False,
#                            keepfiles=False,options=optarg)
#        
#        if outlvl > 0:
#            print(value(time.time() - ts), " s")
#            print(results.solver.message)  
#            print("""    """)   


#%%        # ---------------------------------------------------------------------
        # 5th Initialisation Step - Outlet conditions
        # Initialise variables
        print()
        print("Initialise outlet conditions") 
        
        # Activate constraints
        blk.eq_f7.activate()
        blk.eq_f8.activate()
        blk.eq_f9.activate()
        blk.eq_f10.activate()
        
        # Unfix variables
        blk.Gas_Out_P.unfix()
        blk.Solid_Out_M.unfix()
        blk.Solid_Out_Ts.unfix()
        blk.Solid_Out_x.unfix() 
        
        # Solve for outlet conditions 
        ts = time.time()
        if outlvl > 0:
            print(blk,"- Step 5, initialize oulet conditions terms Time: ", end="")       

        results = opt.solve(blk,tee=True,symbolic_solver_labels=False,
                            keepfiles=False,options=optarg)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)   

'''
