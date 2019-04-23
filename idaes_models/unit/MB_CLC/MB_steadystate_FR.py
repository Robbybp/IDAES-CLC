#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 08:55:31 2019

Moving Bed Reactor Model - based on Kim et al. (2016) and Mondino et al. (2017)
Counter-current (gas phase inlet @ z=0, solid phase inlet @ z=1)
Steady-state

Assumptions: - Ideal Gas behavior
             - Constant solids velocity through the bed
             - Non-isothermal/adiabatic/non-adiabatic/ 
               heat loss through refractory-lined reactor wall
             - Interstitial gas velocity calculated using Ergun Eq. or 
               simplified pressure balance
               
@author: Anca Ostace (aostace)
"""

from __future__ import division
from __future__ import print_function

from pyomo.environ import Param, Var, Constraint, Expression, Set, Reals, \
                            acos, PositiveReals, sqrt, exp, Block, value, \
                            TransformationFactory
from pyomo.dae import DerivativeVar, ContinuousSet

import time
import math
#import importlib

# Import IDAES cores
from idaes_models.core import UnitModel, ProcBlock
import idaes_models.core.util.misc as model_util

from pyomo.opt import SolverFactory

__all__ = ['MB Reactor']

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
        """     
        
        self.dae_method = kwargs.pop("dae_method","OCLR")

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
            
        print(self.fe_set)  
        
#        self.nfe = kwargs.pop("nfe",15)
        self.ncp = kwargs.pop("ncp",3)
        
        
        UnitModel.__init__(self, *args, **kwargs)

    def build(self, *args, **kwargs):
        
        self._make_params_props()
        self._make_params()
        self._make_domain()
        
        self._make_vars()
        self._make_vars_props()
        self._make_derivatives()
        
        self._make_constraints_props()
        self._mass_balance()
        self._energy_balance()
        self._pressure_balance()
        self._make_bdry_conds()
        
        self._transfer_coefficients()
        
        self._dae_transform()
        
        
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
        self.k0 = Param(self.rxn_idx, default = 8e-4,
                        doc = 'pre-exponential factor, mol^(1-N_rxn)m^(3*N_rxn -2)/s')
        self.E = Param(self.rxn_idx, default = 4.9e4,
                       doc = 'Activation energy, kJ/mol')
        self.N_rxn = Param(self.rxn_idx, default = 1.3,
                           doc = 'reaction order in gas species, (-)')
        self.b_rxn = Param(self.rxn_idx, default = 12,
                           doc = 'reaction stoich coeff, (-)')     
        # Available volume for reaction - from EPAT report (1-ep)'
        self.a_vol = Param(within = Reals, default = 0.28, 
                         doc= 'available reaction vol. per vol. of OC')
        #Viscosity constants - Reference: Perry and Green Handbook; McGraw Hill, 2008
        self.mu_param = {('CH4','a'):5.2546e-7,('CH4','b'):0.59006,('CH4','c'):105.67,('CH4','d'):0,
                         ('CO2','a'):2.148e-6,('CO2','b'):0.46,('CO2','c'):290,('CO2','d'):0,
                         ('H2O','a'):1.7096e-8,('H2O','b'):1.1146,('H2O','c'):0,('H2O','d'):0}
        
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
        self.hf_in = Param(default=0.42757,
                        doc='Gas-solid heat transfer coeff. [kJ/(m2Ks)]')
        
        self.hw_in = Param(default=0.0857,
                doc='Internal convective heat transfer coeff., [kJ/(m2Ks)]')
        
        self.Cpw = Param(default=0.5,
                         doc='Heat capacity of steel wall [kJ/kg/K]')
        self.rhow = Param(default=8050,
                doc='Reactor wall denisty [kg/m3]')
        
        # h0, aw, aw1 and Tamb not used when adiabatic
        self.h0 = Param(default=0.28094736,
                doc='Wall-ambiance heat transfer coefficient, [kJ/(m2 K s)]')   
        self.aw = Param(default=3.2129, doc='[1/m]')
        self.aw1 = Param(default=3.2727, doc='[1/m]')
        self.Tamb = Param(default=298.15, doc='Ambient temperature [K]')
        
        self.w_th = Param(default=0.005,
                          doc='Reactor wall thickness [m]')
        self.k_steel = Param(
                default=16.3*1e-3, 
                doc='Thermal conductivity of the steel reactor wall [kW/m/K]')
        self.k_refractory = Param(
                default=1.2*1e-3, 
                doc='Thermal conductivity of the refractory lining [kW/m/K]')
        
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
        self.alfa_air = Param(
                default=22.873433*1e-6,
                doc='Thermal diffusivity of air at the film temp. [m2/s]')       
        
        self.tuning_param = Param(default=1, mutable=True)
        self.tuning_param2 = Param(default=1, mutable=True)
        
        
    def _make_domain(self):
        """
        Create the axial dimensionless domain for the model 
        (along length of reactor)
        """
        self.z = ContinuousSet(bounds=(0.0,1.0), domain=Reals,
                               initialize=self.fe_set,
                               doc='Axial dimensionless domain')      
        
    def _make_vars(self):
        """
        Create the model variables. 
        """
        # Reactor bed dimensions and characteristics
        self.L = Var(
                domain=Reals, bounds=(0.0,100.0), initialize=1.5,
                doc='Reactor bed length [m]')
        self.Dr = Var(
                domain=Reals, bounds=(0.0,40.0), initialize=2,
                doc='Reactor inner diameter [m]')
#        self.w_th = Var(
#                domain=Reals, bounds=(0.0,1.0), initialize=0.005,
#                doc='Reactor wall thickness [m]')
        self.refractory_th = Var(
                domain=Reals, bounds=(0.0,100.0), initialize=0.4,
                doc='Refractory thickness [m]')          
        self.eps = Var(
                domain=Reals, bounds=(0.0,1.0), initialize=0.8,
                doc='Bed voidage [-]')             
#        self.L = Var(
#                domain=Reals, bounds=(0.1,25.0), initialize=1.5,
#                doc='Reactor bed length [m]')
#        self.Dr = Var(
#                domain=Reals, bounds=(0.1,25.0), initialize=2,
#                doc='Reactor inner diameter [m]')
##        self.w_th = Var(
##                domain=Reals, bounds=(0.0,10.0), initialize=0.005,
##                doc='Reactor wall thickness [m]')
#        self.refractory_th = Var(
#                domain=Reals, bounds=(0.0,10.0), initialize=0.4,
#                doc='Reactor wall thickness [m]')          
#        self.eps = Var(
#                domain=Reals, bounds=(0.0,1.0), initialize=0.8,
#                doc='Bed voidage [-]')        
        
        # Gas Inlet Conditions   
#        self.Gas_In_M = Var(
#                domain=Reals, bounds=(0.0,1e5), initialize=100.0,
#                    doc='Inlet gas mass flow rate [kg/s]')
        self.Gas_In_F = Var(
                domain=Reals, bounds=(0.0,1e5), initialize=100.0,
                    doc='Inlet gas mmolar flow rate [mol/s]')
        self.Gas_In_P = Var(
                domain=Reals, bounds=(0.0,1e2), initialize=1.05,
                    doc='Pressure of gas at gas inlet [bar]')
        self.Gas_In_Tg = Var(
                domain=Reals, bounds=(0.0,1000.0), initialize=298.15,
                    doc='Temperature of gas at gas inlet [K]')
        self.Gas_In_y = Var(
                self.GasList, domain=Reals, bounds=(0.0,1.0), initialize=0.1,
                    doc='Mole fractions of gas species at gas inlet [mol/mol]')

        # Gas Outlet Conditions
        self.Gas_Out_P = Var(
                domain=Reals, bounds=(1.0,1e2), initialize=1.05,
                    doc='Pressure of gas at gas outlet [bar]')
        
        # Solid Inlet Conditions   
        self.Solid_In_M = Var(
                domain=Reals, bounds=(0.0,1e5), initialize=2280,
                    doc='Inlet solid mass flow rate [kg/s]')
        self.Solid_In_Ts = Var(
                domain=Reals, bounds=(0.0,2000.0), initialize=293.15,
                    doc='Temperature of gas at gas inlet [K]')
        self.Solid_In_x =  Var(
                self.SolidList, domain=Reals, bounds=(0.0,1.0), initialize=0.1,
                doc='Mass fraction of solid species at solid inlet [kg/kg]')  
        
        # Solid Outlet Conditions   
        self.Solid_Out_M = Var(
                domain=Reals, bounds=(0.0,1e5), initialize=2280,
                    doc='Inlet solid mass flow rate [kg/s]')
        self.Solid_Out_M_Comp = Var(
                self.SolidList, domain=Reals, bounds=(0.0,1e5), 
                initialize=2280,		
                doc='Inlet solid mass flow rate [kg/s]')         
        self.Solid_Out_Ts = Var(
                domain=Reals, bounds=(0.0,2000.0), initialize=293.15,
                    doc='Temperature of gas at gas inlet [K]')
        self.Solid_Out_x =  Var(
                self.SolidList, domain=Reals, bounds=(0.0,1.0), initialize=0.1,
                doc='Mass fraction of solid species at solid inlet [kg/kg]') 
        
        # Gas Phase Concentrations
        self.CgT = Var(
                self.z, domain=Reals, initialize=1.0, bounds=(0.0,1e5),
                doc='Bulk gas total concentration [mol/m3]')
        self.Cg = Var(
                self.z, self.GasList, domain=Reals, initialize=1.0,
                bounds=(0.0,1e5), 
                doc='Bulk gas component concentrations [mol/m3]')
        self.Ctrans = Var(
                self.z, self.GasList, domain=Reals, initialize=0.0, 
                doc='Moles transferred from Gas to Solid [mol/m3/s]')
        self.vgCg = Var(
                self.z, self.GasList, domain=Reals, initialize=1,
                doc="Gas (velocity)*(component concentration) (mol/m2/s)")  
        self.X_gas = Var(
                domain=Reals, bounds=(0.0,1.0), initialize=0.0,
                doc="Methane conversion [-]")

        # Gas Mole Fractions
        self.y = Var(
                self.z, self.GasList, domain=Reals, #bounds=(0.0,1.0),
                initialize=0.1, doc='Gas phase mole fractions [-]')
        self.ytot = Var(
                self.z, domain=Reals, initialize=0.1, #bounds=(0.0,1.0),
                doc='Gas phase total mole fraction [-]')  
        
        # Gas flowrates
        self.Ftotal = Var(
                self.z, domain=Reals, initialize=100.0, bounds=(0.0,1e5),
                doc='Total molar flow rate of gas [mol/s]')
        self.F = Var(
                self.z, self.GasList, domain=Reals, initialize=100.0,
                doc='Molar flow rate of gas components [mol/s]')      
        self.Gas_M = Var(
                self.z, self.GasList, domain=Reals, initialize=100.0,
                doc='Mass flow rate of gas components [kg/s]')
        
        
        # Solid Phase Loading     
        self.qT = Var(
                    self.z, domain=Reals, initialize=1.0, bounds=(0.0,1e5),
                    doc='Total solid loading [kg/m3]')
        self.q = Var(
                    self.z, self.SolidList, domain=Reals, initialize=1.0,
                    bounds=(0.0,1e5), doc='Solid loading [kg/m3]')
        self.qtrans = Var(
                    self.z, self.SolidList, domain=Reals,initialize=0.0,
                    doc='Moles transferred from Gas to Solid [kg/m3/s]')
        self.vsq = Var(
                self.z, self.SolidList, domain=Reals, initialize=1,
                doc="Solids (velocity)*(loading) (kg/m2/s)")  
        self.X_Fe2O3 = Var(
                domain=Reals, bounds=(0.0,1.0), initialize=0.0,
                doc="Oxigen carrier conversion [-]")
        
        # Solid mass fractions
        self.x = Var(
                self.z, self.SolidList, domain=Reals, #bounds=(0.0,1.0),
                initialize=0.1, doc='Solid phase mass fractions [-]')
        self.xtot = Var(
                self.z, domain=Reals, initialize=0.1, #bounds=(0.0,1.0),
                doc='Solid phase total mass fraction [-]')            
                
        # Mass flowrates
        self.Solid_M_total = Var(
                self.z, domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='Total mass flow rate of solids [kg/s]')
        self.Solid_M = Var(
                self.z, self.SolidList, domain=Reals, 
                bounds=(0.0,1e7), initialize=1.0,
                doc='Mass flow rate of solid components [kg/s]') 
        self.Solid_F_total = Var(
                self.z, domain=Reals, bounds=(0.0,1e7), initialize=1.0,
                doc='Total molar flow rate of solids [mol/s]')
        self.Solid_F = Var(
                self.z, self.SolidList, domain=Reals, 
                bounds=(0.0,1e7), initialize=1.0,
                doc='Molar flow rate of solid components [mol/s]') 
        self.mFe_mAl = Var(
                self.z, domain=Reals, initialize=1,
                doc='Fe to Al mass ratio [-]')
        
        # Temperatures and Pressures
        self.P = Var(
                self.z, domain=Reals, bounds=(0.0,1e2), initialize=1.05,
                doc='Pressure [bar]')        
        self.Tg = Var(
                self.z, domain=Reals, bounds=(0.0,2000.0), initialize=298.15,
                doc='Gas phase temperature [K]')
        self.Ts = Var(
                self.z, domain=Reals, bounds=(0.0,2000.0), initialize=298.15,
                doc='Solid phase temperature [K]')           
        self.Tw = Var(
                self.z, domain=Reals, bounds=(0.0,10000.0), initialize=298.15,
                doc='Wall temperature [K]')
        self.Tg_GS = Var(
                self.z, domain=Reals, initialize=0.0)
        self.Tg_GW = Var(
                self.z, domain=Reals, initialize=0.0)  
        self.Tg_refractory = Var(
                self.z, domain=Reals, initialize=0.0)
        self.Ts_dHr= Var(
                self.z, domain=Reals, initialize=0.0)
        self.Tw_GW = Var(
                self.z, domain=Reals, initialize=0.0) 
        self.Tw_Wamb = Var(
                self.z, domain=Reals, initialize=0.0) 
        
        # Velocities
        self.vg = Var(
                  domain=Reals, bounds=(-0.001,1e3), initialize=0.05,
                  doc='Superficial gas velocity [m/s]')  # at inlet
        self.vgz = Var(
                self.z, domain=Reals, bounds=(-0.001,1e3), initialize=0.05,
                doc='Superficial gas velocity [m/s]')
        self.umf = Var(
                self.z, domain=Reals,bounds=(-0.001,1e3), initialize=1.00584,
                doc='Minimum fluidization velocity [m/s]')
        self.vs = Var(
                  domain=Reals, initialize=0.005,
                  doc='Solids velocity (downwards) [m/s]') # not negative  
#        self.vsz = Var(
#                self.z, domain=Reals, initialize=0.005,
#                  doc='Solids velocity (downwards) [m/s]') # not negative        
        self.vel_diff=Var(
                self.z, domain=Reals, bounds=(0.0,1e5),
                doc='Variable to store "umf-vgz", used in optimization')                
        
        # Dimensionless numbers
        self.Rep = Var(
                self.z, domain=Reals, bounds=(0.0,1e8), initialize=1.0,
                doc='Particle Reynolds number (dimensionless)')
        self.Re = Var(
                self.z, domain=Reals, bounds=(0.0,1e8), initialize=1.0,
                doc='Reynolds number (dimensionless)')
        self.Pr = Var(
                self.z, domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='Prandtl number (dimensionless)')
        self.Pr_ext = Var(
                domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc=' External (air) Prandtl number (dimensionless)')
        self.Ra = Var(
                domain=Reals, bounds=(0.0,1e17), initialize=1.0,
                doc=' External (air) Prandtl number (dimensionless)')        
        self.Nu = Var(
                self.z, domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='Nusselt number (dimensionless)')
        self.Nuw = Var(
                self.z, domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='Nusselt number (dimensionless)')
        self.Nu_ext = Var(
                domain=Reals, bounds=(0.0,1e5), initialize=1.0,
                doc='External Nusselt number (dimensionless)')
        self.hf = Var(
                self.z, initialize=0.47502,
                        doc='Gas-solid heat transfer coeff. [kJ/(m2Ks)]')
        self.hw = Var(
                self.z, initialize=0.21343,
                doc='Internal convective heat transfer coeff., [kJ/(m2Ks)]')
        self.hext = Var(
                bounds=(0.0,1e5), initialize=0.0018,
                doc='External natural convective HTC [kJ/(m2Ks)]')
        self.hext2 = Var(
                self.z, bounds=(0.0,1e5), initialize=0.02,
                doc='External natural convective HTC [kJ/(m2Ks)]')
        self.U = Var(
                self.z, bounds=(0.0,1e5), initialize=0.0017,
                doc='Overall heat transfer coeff. (including refractory), [kJ/(m2Ks)]')
        self.Uw = Var(
                self.z, bounds=(0.0,1e5), initialize=0.0017,
                doc='Overall heat transfer coeff., [kJ/(m2Ks)]')
              
        # Dummy variables
#        self.dummy = Var(
#                self.z, domain=Reals, initialize = 1.0, 
#                doc='Dummy variable used to store dP/dz')
    
    
    def _make_vars_props(self):
        """
        Create properties variables. These would belong to a property package, 
        but have been moved here because it somehow improved convergence. Will 
        be moved back out while transitioning to the new framework.
        """ 
#        self.H_comp_s = Var(self.z, self.comp, domain = Reals, initialize = 1.0, 
#                        doc = 'Solid phase component Enthalpy at Ts(K) - J/(mol comp)')
        self.DH_rxn_s = Var(self.z, domain = Reals, initialize = 1.0, 
                        doc = 'Heat of rxn at system T, J/(mol reaction)')
        self.cp_sol = Var(self.z,domain = Reals, initialize = 1.0,
                        doc = 'Heat capacity of solid particles, J/kg.K')
        self.rho_bulk = Var(self.z, domain=Reals, initialize = 1951.22,
                            doc = 'bulk density - kg/m3')
        self.MW_vap = Var(self.z, domain=Reals, initialize = 0.03,
                          doc = 'mol wt. of gas - kg/mol')
        self.rho_vap = Var(self.z, domain=Reals, initialize = 1.0,
                        doc = 'gas mass density - kg/m3')
        self.mu_vap = Var(self.z, domain=Reals, initialize = 1e-5,
                        doc = 'dynamic viscosity of gas, - Pa s = kg/m/s')
        self.cp_gas = Var(self.z, domain=Reals, initialize = 1.0,
                        doc='Gas phase heat capacity [kJ/kg/K]') 
        self.cp_vap = Var(self.z, domain = Reals, initialize = 1.0,
                        doc = 'heat capacity of gas - J/mol/K')
        self.cv_vap = Var(self.z, domain = Reals, initialize = 1.0,
                        doc = 'heat capacity of gas - J/mol/K')
        self.k_cpcv = Var(self.z, domain = Reals, initialize = 1.0,
                        doc = 'specific heat ratio (-)')        
        self.k_vap = Var(self.z, domain = Reals, initialize = 1e-1, 
                        doc = 'thermal conductivity of gas, J/mKs')
        self.X = Var(self.z, domain=Reals, initialize=0.0, 
                        doc = 'fraction of metal oxide converted')
        self.X_term = Var(self.z, domain = Reals, initialize = 1.0, 
                        doc = 'reformulation term for X to help eqn scaling') 
        self.k = Var(self.z, self.rxn_idx, domain=Reals, initialize=1.0, 
                        doc = 'kinetic rate constant')
        self.r_gen = Var(self.z, self.rxn_idx, domain=Reals, initialize = 0.0, 
                        doc = 'gen. rate expression, units= mol_rxn/gOC.s')          
        self.rg = Var(self.z, self.GasList, domain=Reals, initialize = 0.0,
                        doc = 'gcomp. total rate expression, units= mol/m^3.s')        
        self.rs = Var(self.z, self.SolidList, domain=Reals, initialize = 0.0,
                      doc = 'scomp. rate expression, units= kg/m^3.s')
    
    
    def _make_derivatives(self):
        """
        Make derivative variables
        """
        self.dvgCgdz = DerivativeVar(self.vgCg, wrt=self.z)
        self.dvsqdz = DerivativeVar(self.vsq, wrt=self.z)
        
        self.dPdz = DerivativeVar(self.P, wrt=self.z) 
        
        self.dTgdz = DerivativeVar(self.Tg, wrt=self.z)
        self.dTsdz = DerivativeVar(self.Ts, wrt=self.z)        


    def _make_constraints_props(self):
        """ 
        Create property constraints. These would belong to a property package, 
        but have been moved here because it somehow improved convergence. Will 
        be moved back out while transitioning to the new framework.
        """
        def rule_eq_H_comp_s(b,z,i):
            return 1e3*(b.cp_param[i,'a']*(b.Ts[z]/1000)
                        + b.cp_param[i,'b']*((b.Ts[z]/1000)**2)/2
                        + b.cp_param[i,'c']*((b.Ts[z]/1000)**3)/3
                        + b.cp_param[i,'d']*((b.Ts[z]/1000)**4)/4
                        - b.cp_param[i,'e']/(b.Ts[z]/1000) + b.cp_param[i,'f']
                        - b.cp_param[i,'h'])
        self.H_comp_s = Expression(self.z, self.comp, rule=rule_eq_H_comp_s)      

        def rule_eq_h3(b,z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.DH_rxn_s[z] == sum(b.stoic[i]*(b.H_comp_s[z,i] 
                                        + b.Hf_comp[i]) for i in b.comp)
        self.eq_h3 = Constraint(self.z, rule = rule_eq_h3)     
        
        # Solid properties
        def rule_expr_cp_comp_s(b,z,i):
            return b.cp_param[i,'a'] + b.cp_param[i,'b']*(b.Ts[z]/1000) \
                        + b.cp_param[i,'c']*(b.Ts[z]/1000)**2 \
                        + b.cp_param[i,'d']*(b.Ts[z]/1000)**3 \
                        + b.cp_param[i,'e']/((b.Ts[z]/1000)**2)
        self.cp_comp_s = Expression(self.z, self.comp, rule=rule_expr_cp_comp_s) # Solid phase - cp_comp_s units = J/mol.K
        
        def rule_eq_cp_comp_s(b,z):
            return b.cp_sol[z] == sum(b.cp_comp_s[z,j]*(1000/b.MW[j])
                                      *b.x[z,j] for j in self.SolidList)
        self.eq_h1 = Constraint(self.z, rule=rule_eq_cp_comp_s)  

        # Bulk density
        def rule_rho_bulk(b,z):
            return b.rho_bulk[z] == (1-b.eps)*b.rho_sol + b.eps*b.rho_vap[z]
        self.eq_rho_bulk = Constraint(self.z, rule=rule_rho_bulk)
        
        # Gas mixture molecular weight
        def rule_eq_MW_vap(b,z):
            return b.MW_vap[z] == 1e-3*sum(b.y[z,i]*b.MW[i] for i in b.GasList)
        self.eq_q2 = Constraint(self.z, rule=rule_eq_MW_vap)
        # Gas phase density
        def rule_rho_vap(b,z):
            return b.rho_vap[z] == b.MW_vap[z]*b.P[z]/(b.R*1e-5*b.Tg[z])
        self.eq_q3 = Constraint(self.z, rule=rule_rho_vap)
        # Components viscosity 
        def rule_eq_q9(b,z,i):
            return  b.mu_param[i,'a']*(b.Tg[z]**b.mu_param[i,'b']) \
                        / ((1 + (b.mu_param[i,'c']/b.Tg[z])) \
                        + (b.mu_param[i,'d']/(b.Tg[z]**2))) 
        self.mu_comp = Expression(self.z,self.GasList, rule=rule_eq_q9)
        # Viscosity of gas mix
        def rule_eq_q10(b,z):
            return 1e6*b.mu_vap[z] == 1e6*sum(b.y[z,i]*b.mu_comp[z,i] \
                                    /(sum(b.y[z,j]*(b.MW[j]/b.MW[i])**0.5 \
                                    for j in b.GasList)) for i in b.GasList)
        self.eq_q10 = Constraint(self.z, rule=rule_eq_q10)

        # Heat capacity
        def rule_eq_q17(b,z,i):
            return b.cp_param[i,'a'] + b.cp_param[i,'b']*(b.Tg[z]/1000) \
                        + b.cp_param[i,'c']*(b.Tg[z]/1000)**2 \
                        + b.cp_param[i,'d']*(b.Tg[z]/1000)**3 \
                        + b.cp_param[i,'e']/((b.Tg[z]/1000)**2)
        self.cp_comp_g = Expression(self.z, self.GasList, rule=rule_eq_q17)
            
        def rule_eq_q18(b,z):
            return b.cp_vap[z] == sum(b.cp_comp_g[z,j]*b.y[z,j] for j in b.GasList)
        self.eq_q18 = Constraint(self.z, rule=rule_eq_q18)         
        
        def rule_eq_q18a(b,z):
            return b.cp_gas[z]*b.MW_vap[z] == b.cp_vap[z]*1e-3
        self.eq_q18a = Constraint(self.z, rule=rule_eq_q18a)
        
        def rule_eq_q18b(b,z):
            return b.cv_vap[z] == b.cp_vap[z] - b.R
        self.eq_q18b = Constraint(self.z, rule=rule_eq_q18b)

        def rule_eq_q18c(b,z):
            return b.k_cpcv[z] ==  b.cp_vap[z]/b.cv_vap[z]
        self.eq_q18c = Constraint(self.z, rule=rule_eq_q18c)        
                    
        # Thermal conductivity of gas  
        def rule_eq_q6(b,z,i):
            return b.k_param[i,'a']*(b.Tg[z]**b.k_param[i,'b']) \
                        / ((1 + (b.k_param[i,'c']/b.Tg[z])) \
                        + (b.k_param[i,'d']/(b.Tg[z]**2)))                         
        self.k_comp = Expression(self.z, self.GasList, rule=rule_eq_q6)

        def rule_eq_q7(b,z,i,j):
            return (1 + ((b.k_comp[z,i]/b.k_comp[z,i])**0.5) \
                        * ((b.MW[j]/b.MW[i])**0.25))**2 \
                        / (8*(1+(b.MW[j]/b.MW[i])))**0.5    
        self.A_bin = Expression(self.z, self.GasList, self.gcomp, rule=rule_eq_q7)
            
        def rule_eq_q8(b,z):
            return 1e6*self.k_vap[z] == 1e6*sum(b.y[z,i]*b.k_comp[z,i] \
                                    /(sum(b.y[z,j]*b.A_bin[z,i,j]**0.5 \
                                    for j in b.GasList)) for i in b.GasList)
        self.eq_q8 = Constraint(self.z, rule=rule_eq_q8)    
        
        # rate expression constraints - overall rate and species specific rates
        def rule_eq_r1(b, z, i):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.k[z,i] == b.k0[i]*exp(-b.E[i]/(b.R*b.Ts[z]))
        self.eq_r1 = Constraint(self.z, self.rxn_idx, rule=rule_eq_r1, doc = 'kinetic rate constant eqn') 
            
        """
        This equation is for calculating the conv. of fe2 to fe3 in the reactor. It is gotten
        by converting the EPAT conv. equation (based on mole flow) to a mass fraction equivalent.
        
            X =                         x_fe3
                ---------------------------------------------------------
                 x_fe3 + x_fe2*((stoich_fe3/-stoich_fe2)*(MW_fe3/MW_fe2))
                
        where  x_fe2,x_fe3 = wt. fraction of fe2 and fe3 at any
                    given point in the reactor (state variables)                     
        """  

        def rule_eq_r2(b,z):
            return 1e9*b.X[z]*(b.x[z,'Fe3O4'] + (b.MW['Fe3O4']/b.MW['Fe2O3']) \
                        *(b.stoic['Fe3O4']/-b.stoic['Fe2O3'])*b.x[z,'Fe2O3']) \
                        == 1e9*b.x[z,'Fe3O4']
        self.eq_r2 = Constraint(self.z,rule=rule_eq_r2, doc = 'conversion of metal oxide eqn')  

        def rule_eq_r2b(b,z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 1e6*b.X_term[z]**3 == 1e6*(1-b.X[z])**2   
        self.eq_r2b = Constraint(self.z, rule = rule_eq_r2b)        
      
        
        # general reaction rate expression of CH4 with Fe2O3
        def rule_eq_r4(b,z,i):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.r_gen[z,i]*1e6== 1e6*b.scale*b.x[z,'Fe2O3']*(b.a_vol/b.MW['Fe2O3'])\
                                    *3*b.b_rxn[i]*b.k[z,i]*(b.Cg[z,'CH4']**b.N_rxn[i])\
                                    *b.X_term[z]/(b.rhom*b.radg)/(-b.stoic['Fe2O3'])                    
        self.eq_r4 = Constraint(self.z, self.rxn_idx, rule=rule_eq_r4, 
                        doc = 'general rate expression (mole extent basis), \
                                units are mol_reaction/gOC.s')       

        def rule_eq_r5(b,z,i):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.rs[z,i] == b.rho_sol*b.MW[i]*b.stoic[i]*b.r_gen[z,1]
        self.eq_r5 = Constraint(self.z, self.SolidList, rule=rule_eq_r5, 
                                    doc = 'comp specific rate expression, \
                                    units are kgOC/m3/s')

        def rule_eq_r6(b,z,i):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.rg[z,i] == -1e3*b.rho_sol*b.stoic[i]*b.r_gen[z,1]                               
        self.eq_r6 = Constraint(self.z, self.GasList, rule=rule_eq_r6, 
                                    doc = 'comp specific rate expression for a rxn,\
                                    units are mol/m3/s')
            
        
    def _transfer_coefficients(self):
        """
        Calculate dimensionless numbers and heat transfer coefficients.
        """        
        # Particle Reynolds number   #*
        def rule_eq_Rep(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Rep[z]*b.mu_vap[z] == b.vgz[z] \
                                *b.dp*b.rho_vap[z]
        self.eq_Rep = Constraint(self.z, rule=rule_eq_Rep)
        
        # Reynolds number
        def rule_eq_Re(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Re[z]*b.mu_vap[z] == b.vgz[z] \
                                *b.Dr*b.rho_vap[z]
        self.eq_Re = Constraint(self.z, rule=rule_eq_Re)
        
        # Prandtl number - cp_gas [J/kg/K]
        def rule_eq_Pr(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Pr[z] == (b.cp_gas[z]*1e3*
                                       b.mu_vap[z]/
                                       b.k_vap[z])
        self.eq_Pr = Constraint(self.z, rule=rule_eq_Pr)

        # External Prandtl number (for air)
        self.eq_Pr_ext = Constraint(expr=self.Pr_ext*self.k_air == 
                                    self.Cp_air*self.mu_air)
        # External Rayleigh number (for air)
        self.eq_Ra = Constraint(expr=self.Ra*self.nu_air*self.alfa_air == 
                        self.g*self.beta*(self.Tw[1]-self.Tamb)*self.L**3.0)
         # Nusselt number
        def rule_eq_Nu(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Nu[z] == (2.0 + 1.1*(abs(b.Rep[z])**0.6)
                                   *(abs(b.Pr[z])**(1/3))) 
        self.eq_Nu = Constraint(self.z, rule=rule_eq_Nu)
        
        # Nusselt number
        def rule_eq_Nuw(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Nuw[z] == (12.5 + 0.048*b.Re[z])
        self.eq_Nuw = Constraint(self.z, rule=rule_eq_Nuw)

        # Exterior natural convection Nusselt number
        self.eq_Nu_ext = Constraint(
                expr=self.Nu_ext == 0.68 + (0.67*self.Ra**(1/4))
                /((1+(0.492/abs(self.Pr_ext))**(9/12))**(4/9)))
        
        # Gas-solid heat transfer coefficient
        def rule_eq_hf(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 1e6*b.hf[z]*b.dp == \
                                    1e6*b.Nu[z]*b.k_vap[z]*1e-3
        self.eq_hf = Constraint(self.z, rule=rule_eq_hf)
        
        # Gas-wall heat transfer coefficient
        def rule_eq_hw(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 1e6*b.hw[z]*b.eps*b.Dr == \
                                    1e6*b.Nuw[z]*b.k_vap[z]*1e-3
        self.eq_hw = Constraint(self.z, rule=rule_eq_hw)

        # Exterior natural convection Nusselt number
        self.eq_hext = Constraint(expr=self.hext*self.L == 
                                  self.k_air*self.Nu_ext)
        
        # Exterior natural convection Nusselt number
        def rule_eq_hext2(b,z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.hext2[z]*b.L/b.nfe == b.k_air*b.Nu_ext    
        self.eq_hext2 = Constraint(self.z, rule=rule_eq_hext2)
        
        # Overall heat transfer coefficient (to use w/ refractory lining)
        def rule_eq_U(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 1e8*b.U[z] == 1e8/((1/b.hw[z]) \
                            + b.refractory_th/b.k_refractory \
                            + b.w_th/b.k_steel \
                            + (1/b.hext2[z]))
        self.eq_U = Constraint(self.z, rule=rule_eq_U)

        # Overall heat transfer coefficient (wall-ambiance)
        def rule_eq_Uw(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 1e8*b.Uw[z] == 1e8/((1/b.hw[z]) \
                            + b.w_th/b.k_steel \
                            + (1/b.hext2[z]))
        self.eq_Uw = Constraint(self.z, rule=rule_eq_Uw)
                
        
    #========================================================================== 
    def _mass_balance(self):
        """
        Add the mass balance constraints. The boundary conditions are defined 
        under '_make_bdry_conds'
        """
        
        # Bulk gas component mass balance
        def rule_eq_vgCg(b, z, j):    # j indexes GasList, z indexes dz
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 0 == -b.dvgCgdz[z,j] - \
                            (1-b.eps)*b.Ctrans[z,j]*b.L 
        self.eq_vgCg = Constraint(self.z, self.GasList, rule=rule_eq_vgCg)
        
        def rule_eq_vsq(b, z, j): 
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 0 == b.dvsqdz[z,j] + (1-b.eps)*b.qtrans[z,j]*b.L               
        self.eq_vsq = Constraint(self.z, self.SolidList, rule=rule_eq_vsq)
        
        # Moles of gas reacting (transferred from Gas to Solid) [mol/m3/s]
        def rule_eq_Ctransferred(b, z, j):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 1e3*b.Ctrans[z,j] == 1e3*b.rg[z,j]
        self.eq_Ctrans = Constraint(self.z, self.GasList, 
                                        rule=rule_eq_Ctransferred)      
        
        # Moles of solid reacting (transferred from Gas to Solid) [kg/m3/s]
        def rule_eq_qtransferred(b, z, j):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return 1e3*b.qtrans[z,j] == 1e3*b.rs[z,j]
        self.eq_qtrans = Constraint(self.z, self.SolidList, 
                                        rule=rule_eq_qtransferred)
    
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Calculate the total and component GAS molar flow rates
        def rule_eq_F(b,z,j):
            return b.F[z,j] == (b.pi*(b.Dr/2)**2)*b.vgCg[z,j]
        self.eq_F = Constraint(self.z, self.GasList, rule=rule_eq_F)
        
        def rule_eq_Gas_M(b,z,j):
            return b.Gas_M[z,j] == b.F[z,j]*b.MW[j]*1e-3
        self.eq_Gas_M = Constraint(self.z, self.GasList, rule=rule_eq_Gas_M)
        
        def rule_eq_Ftotal(b,z):
            return b.Ftotal[z] == sum(b.F[z,j] for j in b.GasList)
        self.eq_Ftotal = Constraint(self.z, rule=rule_eq_Ftotal)

        # Calculate the gas components concetrations
        def rule_eq_Cg(b,z,j):
            if z == b.z.first():
                return b.Cg[z,j]*b.R*1e-5*b.Gas_In_Tg == \
                                            b.Gas_In_y[j]*b.Gas_In_P
            else:
                return b.Cg[z,j] == b.F[z,j]/((b.pi*(b.Dr/2)**2)*b.vgz[z])  
        self.eq_Cg = Constraint(self.z, self.GasList, rule=rule_eq_Cg)        
                
        # Bulk gas total concentration
        def rule_eq_CgT(b, z):
            if z == b.z.first():
                return b.CgT[z]*b.R*1e-5*b.Gas_In_Tg == b.Gas_In_P 
            else:
                return b.CgT[z] == b.Ftotal[z]/((b.pi*(b.Dr/2)**2)*b.vgz[z])
        self.eq_CgT = Constraint(self.z, rule=rule_eq_CgT)
        
        # Calculate the gas component mole fractions 
        def rule_eq_y(b, z, j):
            return b.y[z,j]*b.CgT[z] == b.Cg[z,j]
        self.eq_y = Constraint(self.z, self.GasList, rule=rule_eq_y)
        
        def rule_eq_ytot(b,z):
            return b.ytot[z] == sum(b.y[z,j] for j in b.GasList) 
        self.eq_ytot = Constraint(self.z, rule=rule_eq_ytot) 
        
        self.eq_gas_conversion = Constraint(expr=self.X_gas == 
                                    (1 - self.Gas_M[1,'CH4']/self.Gas_M[0,'CH4']))
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Calculate the total and component SOLID mass flow rates
        def rule_eq_Solid_M(b,z,j):
            return b.Solid_M[z,j] == (b.pi*(b.Dr/2)**2)*b.vsq[z,j]
        self.eq_Solid_M = Constraint(self.z,self.SolidList,rule=rule_eq_Solid_M)
        
        def rule_eq_Solid_M_total(b,z):
            return b.Solid_M_total[z] == sum(b.Solid_M[z,j] \
                                              for j in b.SolidList)
        self.eq_Solid_M_total = Constraint(self.z, rule=rule_eq_Solid_M_total)
        
        def rule_eq_Solid_F(b,z,j):
            return b.Solid_F[z,j] == b.Solid_M[z,j]/(b.MW[j]*1e-3)
        self.eq_Solid_F = Constraint(self.z,self.SolidList,rule=rule_eq_Solid_F)
        
        def rule_eq_Solid_F_total(b,z):
            return b.Solid_F_total[z] == sum(b.Solid_F[z,j] \
                                              for j in b.SolidList)
        self.eq_Solid_F_total = Constraint(self.z, rule=rule_eq_Solid_F_total)
        
        # Total solid loading [kg/m3]
        def rule_eq_qT(b, z):
            return b.qT[z] == sum(b.vsq[z,j]/b.vs for j in b.SolidList) 
        self.eq_qT = Constraint(self.z, rule=rule_eq_qT)
        
        # Calculate the solid components loading
        def rule_eq_q(b,z,j):
            return b.q[z,j] == b.vsq[z,j]/b.vs 
        self.eq_q = Constraint(self.z, self.SolidList, rule=rule_eq_q)        
        
        # Calculate the solid phase mass fractions 
        def rule_eq_x(b, z, j):
            return b.x[z,j]*b.qT[z] == b.q[z,j]  
        self.eq_x = Constraint(self.z, self.SolidList, rule=rule_eq_x)

        def rule_OC_spec(b,z):
            return ((b.Solid_M[z,'Fe2O3']+b.Solid_M[z,'Fe3O4'])*0.6994
                                /(b.Solid_M[z,'Al2O3']*0.529227) 
                                == b.mFe_mAl[z])
        self.OC_spec = Constraint(self.z,rule=rule_OC_spec)
        
        def rule_eq_xtot(b,z):
            return b.xtot[z] == sum(b.x[z,j] for j in b.SolidList) 
        self.eq_xtot = Constraint(self.z, rule=rule_eq_xtot)  
        
        self.eq_Fe2O3_conversion = Constraint(expr=self.X_Fe2O3 == 
                        1 - self.Solid_M[0,'Fe2O3']/self.Solid_M[1,'Fe2O3'])
        
        
    #==========================================================================     
    def _energy_balance(self):
        """ 
        Add the energy balance constraints. 
        """
        
        def rule_eq_Tg(b, z):
            if z == b.z.first():
                return b.Tg[z] == b.Gas_In_Tg  #Constraint.Skip 
            else:
                return 0 == - b.rho_vap[z]*b.cp_gas[z]*b.vgz[z] \
                            *b.dTgdz[z] \
                            - b.Tg_GS[z]*b.L - b.Tg_GW[z]*b.L \
                            - b.Tg_refractory[z]*b.L
        self.eq_Tg = Constraint(self.z, rule=rule_eq_Tg)    
        
        def rule_eq_Tg_GS(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Tg_GS[z] == b.tuning_param2 \
                              *(1-b.eps)*b.hf[z]*(b.Tg[z]-b.Ts[z]) \
                              *6/b.dp
        self.eq_Tg_GS = Constraint(self.z, rule=rule_eq_Tg_GS)
        
        def rule_eq_Tg_GW(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Tg_GW[z] == b.hw[z]*(b.Tg[z] - b.Tw[z])*4/b.Dr
        self.eq_Tg_GW = Constraint(self.z, rule=rule_eq_Tg_GW)        

        def rule_eq_Tg_refractory(b, z):
            return b.Tg_refractory[z] == b.U[z]*(b.Tg[z] - b.Tw[z])*4/b.Dr
        self.eq_Tg_refractory = Constraint(self.z, rule=rule_eq_Tg_refractory)        


        def rule_eq_Ts(b, z):
            if z == b.z.first():
                return Constraint.Skip# The BC for Ts is under '_make_bdry_conds' 
            else:
                return 0 == b.rho_sol \
                            *b.cp_sol[z]*1e-3 \
                            *b.vs*b.dTsdz[z] \
                            + b.Tg_GS[z]*b.L \
                            + b.Ts_dHr[z]*b.L
        self.eq_Ts = Constraint(self.z, rule=rule_eq_Ts)    
               
        def rule_eq_Ts_dHr(b, z):
            if z == b.z.first():
                return Constraint.Skip 
            else:
                return b.Ts_dHr[z] == b.tuning_param \
                                  *(1-b.eps)*b.rho_sol \
                                  *(-b.DH_rxn_s[z]) \
                                  *b.r_gen[z,1]
        self.eq_Ts_dHr = Constraint(self.z, rule=rule_eq_Ts_dHr)
        
        def rule_eq_Tw(b, z):
            return 0 == b.Tw_GW[z] - b.Tw_Wamb[z]
        self.eq_Tw = Constraint(self.z, rule=rule_eq_Tw)
        
        def rule_eq_Tw_GW(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Tw_GW[z] == b.aw*b.hw[z]*(b.Tg[z]-b.Tw[z])
        self.eq_Tw_GW = Constraint(self.z, rule=rule_eq_Tw_GW)
        
        def rule_eq_Tw_Wamb(b, z):
            if z == b.z.first():
                return Constraint.Skip
            else:
                return b.Tw_Wamb[z] == b.aw1*b.Uw[z]*(b.Tw[z]-b.Tamb)
            #/(b.rhow*b.Cpw)
        self.eq_Tw_Wamb = Constraint(self.z, rule=rule_eq_Tw_Wamb)
        
        
    #==========================================================================  
    def _pressure_balance(self):
        """ 
        Add the pressure balance constraints, Ergun equation.
        """
        
        # Ideal Gas Law
        def rule_eq_P(b,z):
            if z == 0:
                return b.P[z] == b.Gas_In_P  #+ b.ErgunRHS
            else:
                return b.P[z] == b.CgT[z]*(b.R*1e-5)*b.Tg[z]
        self.eq_P = Constraint(self.z, rule=rule_eq_P)
        
        # Simplified pressure balance used to calculate gas velocity - Mondino et al. (2017)
        def rule_simplifiedP(b,z):
            if z == b.z.first():
                return b.vgz[z] == b.vg     # Inlet velocity
            else:
                return -b.dPdz[z]*1e5 == (
                                (b.rho_sol - b.rho_vap[z])
                                *0.2*b.vgz[z]*b.L)
        self.eq_simplifiedP = Constraint(self.z, rule=rule_simplifiedP)

        # Ergun equation used to calculate gas velocity - does not always converge
        def rule_gasvelocity(b,z):
            if z == b.z.first():
                return b.vgz[z] == b.vg     # Inlet velocity
            else:
                return -b.dPdz[z]*1e5 == ((150.0*b.mu_vap[z] 
                                *(1-b.eps)**2*(b.vgz[z] + b.vs) 
                                /(b.dp**2*b.eps**3)) 
                                + (1.75*b.rho_vap[z]*(1-b.eps) 
                                *abs(b.vgz[z] + b.vs)*(b.vgz[z] + b.vs) 
                                /(b.dp*b.eps**3)))*b.L
        self.eq_ErgunP = Constraint(self.z, rule=rule_gasvelocity)
        
        # Compute the minimum fluidized velocity 
        def rule_minfluidvel(b,z):
            return (1.75/(0.45)**3)*(b.dp*b.umf[z] 
                       *b.rho_vap[z]/b.mu_vap[z])**2 \
                       + (150*(1-0.45)/(0.45)**3)*(b.dp \
                       *b.umf[z]*b.rho_vap[z]/b.mu_vap[z]) == \
                       b.dp**3*b.rho_vap[z] \
                       *(b.rho_sol-b.rho_vap[z])*b.g \
                       /(b.mu_vap[z]**2)
        self.eq_minfluidvel = Constraint(self.z, rule=rule_minfluidvel)               

        # Store the difference between the minimum fluidization velocity and 
        # gas velocity, to use as either post-solve check or in optimization
        def rule_velocity(b,z):
            return b.vel_diff[z] == b.umf[z] - b.vgz[z]
        self.velocity_diff = Constraint(self.z, rule=rule_velocity)
        
        
    #==========================================================================
    def _make_bdry_conds(self):
        """
        Boundary conditions for balance equations.
        And inlet velocity of gas and solids.
        """
        
        # BC for gas components mass balance
        def rule_eq_vgCg_BC(b,j):
            return 1e2*b.vgCg[0,j] == 1e2*b.eps*b.Gas_In_F*b.Gas_In_y[j]\
                                            /(b.eps*b.pi*(b.Dr/2.0)**2*1)
        self.eq_vgCg_BC = Constraint(self.GasList, rule=rule_eq_vgCg_BC)
        
        # BC for solid components mass balance
        def rule_eq_vsq_BC(b,j):
            return 1e2*b.vsq[1,j] == 1e2*b.Solid_In_M*b.Solid_In_x[j] \
                                            /(b.pi*(b.Dr/2.0)**2*1)
#             Vbed = (1-eps)*Abed*L, but L=1 here because of scaling
        self.eq_vsq_BC = Constraint(self.SolidList, rule=rule_eq_vsq_BC)  
        
        # BC for solid phase energy balance
        self.eq_Ts_BC = Constraint(expr=self.Ts[1] == self.Solid_In_Ts)

        # ---------------------------------------------------------------------
        # Some inlet conditions

        # Inlet velocities    
        self.eq_vg = Constraint(expr=self.vg == self.eps*self.Gas_In_F 
                        *self.MW_vap[0] 
                        /(self.eps*(self.pi*(self.Dr/2.0)**2)*self.rho_vap[0]))  
        self.eq_vs = Constraint(expr=self.vs == self.Solid_In_M 
                            /((self.pi*(self.Dr/2.0)**2)*self.rho_bulk[0]))      

        # ---------------------------------------------------------------------
        # Gas outlet conditions - used when connecting to other unit models
        self.eq_Gas_Out_P = Constraint(expr=self.Gas_Out_P == self.P[1])
        
        # Solids outlet conditions - used when connecting to other unit models
        self.eq_Solid_Out_M = Constraint(expr= self.Solid_Out_M 
                                         == self.Solid_M_total[0])    
        self.eq_Solid_Out_Ts = Constraint(expr= self.Solid_Out_Ts 
                                         == self.Ts[0]) 
        def rule_Solid_Out_x(b,j):
            return b.Solid_Out_x[j] == b.x[0,j]
        self.eq_Solid_Out_x = Constraint(self.SolidList, rule=rule_Solid_Out_x)
        

    def _dae_transform(self):
        """
        Apply DAE transformation to domain.
        """
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
            
            
    def _initialize(blk, y=None, outlvl=0, optarg=None):
        ''' Initialisation routine for unit (default solver ipopt)'''
        
        # Create a solver
        stee = True
        opt = SolverFactory('ipopt')
        opt.options = {'tol': 1e-8,
                       'linear_solver'  : 'ma27',
#                       'mu_init': 1e-8,
#                       'bound_push': 1e-8,
                       'max_cpu_time': 600,
                       'print_level': 5,}
        
        # Strip bounds
        blk.strip_bounds()
        
        # Set solver options
        if outlvl > 1:
            stee = True
        else:
            stee = False

        if optarg == None:
            sopt = {"tol"            : 1e-8,
                    "max_cpu_time"   : 600,
                    "print_level"    : 5}
        else:
            sopt=optarg        
            
        if outlvl > 0:
            print("\n")
            print("----------------------------------------------------------")
            print("MB Fuel Reactor Initialization ...")
            print("\n")
        
        # ---------------------------------------------------------------------
        # Step 1 - Solve advection only in mass balances only  

        # Fix the temperatures
        for z in blk.z:
            blk.Tg[z].fix(value(blk.Gas_In_Tg))
            blk.Ts[z].fix(value(blk.Solid_In_Ts))    
            blk.Tw[z].fix(value(790)) #blk.Gas_In_Tg))
#            blk.Tw[z].fix(373.15)   # External wall temp is wanted to be 100C
        blk.eq_Tg.deactivate()
        blk.eq_Ts.deactivate()
        blk.eq_Tw.deactivate()
        blk.eq_Ts_BC.deactivate()

        for z in blk.z:
            blk.Tg_GS[z].fix(0.0)
            blk.Tg_GW[z].fix(0.0)
            blk.Tg_refractory.fix(0.0)
            blk.Ts_dHr[z].fix(0.0)  
            blk.Tw_GW.fix(0.0)
            blk.Tw_Wamb.fix(0.0)
        blk.eq_Tg_GS.deactivate()
        blk.eq_Tg_GW.deactivate()
        blk.eq_Tg_refractory.deactivate()
        blk.eq_Ts_dHr.deactivate()
        blk.eq_Tw_GW.deactivate()
        blk.eq_Tw_Wamb.deactivate()
        
        for z in blk.z:
            blk.hf[z].fix(0.5)
            blk.hw[z].fix(0.08)
            blk.U[z].fix(0.0017)
            blk.Uw[z].fix(0.0017)
        blk.eq_hf.deactivate()
        blk.eq_hw.deactivate()
        blk.eq_U.deactivate()
        blk.eq_Uw.deactivate()
        
        # Fix the pressure
        for z in blk.z:
            blk.P[z].fix(value(blk.Gas_In_P))   
        blk.eq_P.deactivate()
        
        # Fix the velocity  
        for z in blk.z:
            blk.vgz[z].fix(value(blk.Gas_In_F)*value(blk.MW_vap[0])
                   /(value(blk.eps)*(value(blk.pi)*(value(blk.Dr)/2)**2)
                   *value(blk.rho_vap[0])))    
        blk.eq_simplifiedP.deactivate()     
        blk.eq_ErgunP.deactivate()
        
        # Fix the minimum fluidization velocity
        for z in blk.z:
            blk.umf[z].fix(1.3801)
        blk.eq_minfluidvel.deactivate()
        
        # Fix flow rates
        for z in blk.z:
            for j in blk.GasList:
                blk.F[z,j].fix(value(blk.Gas_In_F)*value(blk.Gas_In_y[j]))
                blk.Gas_M[z,j].fix(value(blk.Gas_In_F)*value(blk.MW_vap[0])
                                                *value(blk.Gas_In_y[j]))
        blk.eq_F.deactivate()
        blk.eq_Gas_M.deactivate()
        
        for z in blk.z:
            blk.Ftotal[z].fix(value(blk.Gas_In_F))
        blk.eq_Ftotal.deactivate()
        
        for z in blk.z:
            for j in blk.SolidList:
                blk.Solid_M[z,j].fix(value(blk.Solid_In_M)
                                            *value(blk.Solid_In_x[j]))
                blk.Solid_F[z,j].fix(value(blk.Solid_M[z,j]))
        blk.eq_Solid_M.deactivate()
        blk.eq_Solid_F.deactivate()
        
        for z in blk.z:
            blk.Solid_M_total[z].fix(value(blk.Solid_In_M))
            blk.Solid_F_total[z].fix(value(blk.Solid_In_M))
        blk.eq_Solid_M_total.deactivate()
        blk.eq_Solid_F_total.deactivate()
        
        # Fix gas phase concentrations
        for z in blk.z:
            for j in blk.GasList:
                blk.Cg[z,j].fix(value(blk.rho_vap[0])*value(blk.Gas_In_y[j]) \
                              /value(blk.MW_vap[0])) 
        blk.eq_Cg.deactivate()
                
        for z in blk.z:
            blk.CgT[z].fix(value(blk.rho_vap[0])/value(blk.MW_vap[0])) 
        blk.eq_CgT.deactivate()
        
        # Fix solid phase concentrations
        for z in blk.z:
            for j in blk.SolidList:
                blk.q[z,j].fix(value(blk.rho_sol)*value(blk.Solid_In_x[j]))
        blk.eq_q.deactivate()
        
        for z in blk.z:
            blk.qT[z].fix(value(blk.rho_sol)) 
        blk.eq_qT.deactivate()
        
        # Fix mole fractions
        for z in blk.z:
            for j in blk.GasList:
                blk.y[z,j].fix(value(blk.Gas_In_y[j]))
        blk.eq_y.deactivate()
        
        for z in blk.z:
            blk.ytot[z].fix(1)
        blk.eq_ytot.deactivate()
        
        # Fix mass fractions
        for z in blk.z:
            for j in blk.SolidList:
                blk.x[z,j].fix(value(blk.Solid_In_x[j]))
        blk.eq_x.deactivate()
        
        for z in blk.z:
            blk.xtot[z].fix(1)
        blk.eq_xtot.deactivate()       
        
        blk.X_gas.fix(0.0)
        blk.eq_gas_conversion.deactivate()
        
        blk.X_Fe2O3.fix(0.0)
        blk.eq_Fe2O3_conversion.deactivate()
        
        # Fix number of transferred moles/kg and deactivate constraint               
        for z in blk.z:
            for j in blk.GasList:
                blk.Ctrans[z,j].fix(0.0)
            for i in blk.SolidList:
                blk.qtrans[z,i].fix(0.0)    
        blk.eq_Ctrans.deactivate()
        blk.eq_qtrans.deactivate()         
        
        # Set initial values of state variables 
        for z in blk.z:
            for j in blk.GasList:
                blk.vgCg[z,j] = value(blk.Gas_In_F)*value(blk.Gas_In_y[j]) \
                                    /(value(blk.pi)*(value(blk.Dr)/2)**2)
                                    
        for z in blk.z:
            for j in blk.SolidList:
                blk.vsq[z,j] = value(blk.Solid_In_M)*value(blk.Solid_In_x[j]) \
                               /(value(blk.pi)*((value(blk.Dr)/2.0)**2)*1) 
                               

#        for z in blk.z:
#            blk.Rep.fix(10)  
#            blk.Re[z].fix(40000)
#            blk.Pr[z].fix(7)
#            blk.Nu[z].fix(4.8)
#            blk.Nuw[z].fix(2110)
#        blk.eq_Rep.deactivate()         
#        blk.eq_Re.deactivate() 
#        blk.eq_Pr.deactivate()
#        blk.eq_Nu.deactivate()
#        blk.eq_Nuw.deactivate()

#        for z in blk.z:
#            blk.dummy[z].fix(0.0)
#        blk.eq_dummy.deactivate()                      


        # Solve mass balance equations with transfer off 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 1, Initialize advection, Time: ", end="")       

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)        

        # ---------------------------------------------------------------------         
        # Step 2 - Turn on flow rates and concentration, mole fraction 
        # calculations

        # Restore bounds
        blk.restore_bounds()
        
        # Unfix gas flow rates
        blk.F.unfix()
        blk.eq_F.activate()
        
        blk.Gas_M.unfix()
        blk.eq_Gas_M.activate()
        
        blk.Ftotal.unfix()
        blk.eq_Ftotal.activate()
        
        blk.Solid_M.unfix()
        blk.eq_Solid_M.activate()
        
        blk.Solid_M_total.unfix()
        blk.eq_Solid_M_total.activate()
        
        blk.Solid_F.unfix()
        blk.eq_Solid_F.activate()
        
        blk.Solid_F_total.unfix()
        blk.eq_Solid_F_total.activate()
        
        # Solve flow rates - transfer off 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 2, Solve flow rates, Time: ", end="")        
            
        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)                  
        
        # ---------------------------------------------------------------------
        # Step 3 - Solve concentrations
        
        # Unfix gas phase concentrations
        blk.Cg.unfix()
        blk.eq_Cg.activate()
        
        blk.CgT.unfix()
        blk.eq_CgT.activate()
        
        # Unfix solid phase concentrations
        blk.q.unfix()
        blk.eq_q.activate()
        
        blk.qT.unfix()
        blk.eq_qT.activate()
                
        # Solve for concentrations - transfer off 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 3, Solve concentrations, Time: ", end="")         

        results = opt.solve(blk, symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)             
        
        # ---------------------------------------------------------------------
        # Step 4 - Solve molar and mass fractions
        
        # Unfix molar fractions
        blk.y.unfix()
        blk.eq_y.activate()
        
        blk.ytot.unfix()
        blk.eq_ytot.activate()
        
        # Unfix mass fractions
        blk.x.unfix()
        blk.eq_x.activate()
        
        blk.xtot.unfix()
        blk.eq_xtot.activate()        
        
                                      
        # Solve mass balance equations with transfer off  
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 4, Calculate molar and mass fractions, Time: ", end="")         

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)           
          
        # ---------------------------------------------------------------------
        # Step 5 - Solve pressures 
        
        # Activate pressure balances
        blk.P.unfix()
        blk.eq_P.activate()
        
        # Unfix and compute the velocity
        blk.vgz.unfix()   
        blk.eq_simplifiedP.activate()         
#        blk.eq_ErgunP.activate()  

        # Solve mass balance equations with transfer off  
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 5, Calculate pressure and initilize velocity, Time: ", end="")         

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)  
        
        # ---------------------------------------------------------------------
#        # Step 5a - Solve gas velocity
#        # Unfix and compute the velocity
#        
#        blk.eq_simplifiedP.deactivate()            
#        blk.eq_ErgunP.activate()  
#
#        # Solve mass balance equations with transfer off  
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,"- Step 5a, Calculate gas velocity using Ergun...")        
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(results.solver.message)  
#            print("""    """)  
            
        # ---------------------------------------------------------------------
        # Step 5b - Solve for minimum fluidization velocity
        
        # Unfix umf
        blk.umf.unfix()
        blk.eq_minfluidvel.activate()
 
        # Solve mass balance equations with transfer off 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 5b, Calculate the minimum fluidization velocity, Time: ", end="")         

        results = opt.solve(blk, symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)  
            
        # ---------------------------------------------------------------------
        # Step 6 - Turn on mass transfer w/ 1% reaction             
        
        # Unfix number of transferred moles/kg and activate constraint  
        blk.Ctrans.unfix()
        blk.eq_Ctrans.activate()
                     
        blk.qtrans.unfix()  
        blk.eq_qtrans.activate()      
        
        for z in blk.z:
            blk.scale = 0.01
        
        # Solve mass balance equations with transfer on 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 6, Turn on mass transfer w/ 1% reaction, Time: ", end="")        

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)  

        # ---------------------------------------------------------------------
        # Step 6a - Turn on mass transfer w/ 10% reaction                  
        
        for z in blk.z:
            blk.scale = 0.1
        
        # Solve mass balance equations with transfer on 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 6a, Turn on mass transfer w/ 10% reaction, Time: ", end="")         

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)  

        # ---------------------------------------------------------------------
        # Step 6b - Turn on mass transfer w/ 20% reaction                 
        
        for z in blk.z:
            blk.scale = 0.2
        
        # Solve mass balance equations with transfer on 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 6b, Turn on mass transfer w/ 20% reaction, Time: ", end="")         

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)  
            
#        # ---------------------------------------------------------------------
#        # Step 6c - Turn on mass transfer w/ 40% reaction                 
#        
#        for z in blk.z:
#            blk.scale = 0.4
#        
#        # Solve mass balance equations with transfer on 
#        ts = time.time()
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,"- Step 6c, Turn on mass transfer w/ 40% reaction, Time: ", end="")         
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(value(time.time() - ts), " s")
#            print(results.solver.message)  
#            print("""    """)  
#
#        # ---------------------------------------------------------------------
#        # Step 6d - Turn on mass transfer w/ 50% reaction                 
#        
#        for z in blk.z:
#            blk.scale = 0.5
#        
#        # Solve mass balance equations with transfer on 
#        ts = time.time()
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,"- Step 6e, Turn on mass transfer w/ 50% reaction, Time: ", end="")         
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(value(time.time() - ts), " s")
#            print(results.solver.message)  
#            print("""    """)  
#                        
#        # ---------------------------------------------------------------------
#        # Step 6e - Turn on mass transfer w/ 60% reaction             
#
#        for z in blk.z:
#            blk.scale = 0.6
#        
#        # Solve mass balance equations with transfer on 
#        ts = time.time()
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,"- Step 6e, Turn on mass transfer w/ 60% reaction, Time: ", end="")         
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(value(time.time() - ts), " s")
#            print(results.solver.message)  
#            print("""    """)  
#
#        # ---------------------------------------------------------------------
#        # Step 6f - Turn on mass transfer w/ 70% reaction             
#
#        for z in blk.z:
#            blk.scale = 0.7
#        
#        # Solve mass balance equations with transfer on 
#        ts = time.time()
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,"- Step 6f, Turn on mass transfer w/ 70% reaction, Time: ", end="")         
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(value(time.time() - ts), " s")
#            print(results.solver.message)  
#            print("""    """)  
# 
#        # ---------------------------------------------------------------------
#        # Step 6g - Turn on mass transfer w/ 80% reaction             
#
#        for z in blk.z:
#            blk.scale = 0.8
#        
#        # Solve mass balance equations with transfer on 
#        ts = time.time()
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,"- Step 6g, Turn on mass transfer w/ 80% reaction, Time: ", end="")         
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(value(time.time() - ts), " s")
#            print(results.solver.message)  
#            print("""    """)  
#            
#        # ---------------------------------------------------------------------
#        # Step 6h - Turn on mass transfer w/ 90% reaction             
#
#        for z in blk.z:
#            blk.scale = 0.9
#        
#        # Solve mass balance equations with transfer on 
#        ts = time.time()
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,"- Step 6h, Turn on mass transfer w/ 90% reaction, Time: ", end="")         
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(value(time.time() - ts), " s")
#            print(results.solver.message)  
#            print("""    """)  
#                                   
        # ---------------------------------------------------------------------
        # Step 6i - Turn on mass transfer w/ full reaction             

        for z in blk.z:
            blk.scale = 1
            
        # Compute conversions    
        blk.X_gas.unfix()
        blk.eq_gas_conversion.activate()
        
        blk.X_Fe2O3.unfix()
        blk.eq_Fe2O3_conversion.activate()            
        
        # Solve mass balance equations with transfer on 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 6i, Turn on mass transfer w/ full reaction, Time: ", end="")         

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)  
             
        # ---------------------------------------------------------------------
        # Step 7 - Turn on heat advection              
        
        # Unfix gas and solids temperature and activate constraint  

        blk.Tg.unfix()
        blk.eq_Tg.activate()
 
        blk.Ts.unfix()
        blk.eq_Ts.activate()
        blk.eq_Ts_BC.activate()    
        
        # Solve heat advection 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 7, Turn on heat convection, Time: ", end="")        
            
        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)  


        # ---------------------------------------------------------------------
        # Step 8 - Turn on 1% reaction heat            
        
        # Set value of tuning parameter, to turn down reaction heat
        blk.tuning_param = 0.01
        
        # Unfix reaction heat term and activate constraint  
        blk.Ts_dHr.unfix()
        blk.eq_Ts_dHr.activate() 
        
        # Solve complete energy balance 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 8, Turn on 1% reaction heat, Time: ", end="")       

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)            

        # ---------------------------------------------------------------------
        # Step 8a - Turn on 10% reaction heat            
        
        # Set value of tuning parameter, to turn down reaction heat
        blk.tuning_param = 0.1
        
        # Unfix reaction heat term and activate constraint  
        blk.Ts_dHr.unfix()
        blk.eq_Ts_dHr.activate() 
        
        # Solve complete energy balance 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 8a, Turn on 10% reaction heat, Time: ", end="")         

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)  
# 
        # ---------------------------------------------------------------------
        # Step 8b - Turn on 25% reaction heat            
        
        # Set value of tuning parameter, to turn down reaction heat
        blk.tuning_param = 0.25
        
        # Unfix reaction heat term and activate constraint  
        blk.Ts_dHr.unfix()
        blk.eq_Ts_dHr.activate() 
        
        # Solve complete energy balance 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 8b, Turn on 25% reaction heat, Time: ", end="")        

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)    
#
#        # ---------------------------------------------------------------------
#        # Step 8c - Turn on 50% reaction heat            
#        
#        # Set value of tuning parameter, to turn down reaction heat
#        blk.tuning_param = 0.5
#        
#        # Unfix reaction heat term and activate constraint  
#        blk.Ts_dHr.unfix()
#        blk.eq_Ts_dHr.activate() 
#        
#        # Solve complete energy balance 
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,"- Step 8c, Turn on 50% reaction heat ...")        
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(results.solver.message)  
#            print("""    """)   
            
        # ---------------------------------------------------------------------
        # Step 8d - Turn on reaction heat      
        
        # Set value of tuning parameter, to turn up reaction heat
        blk.tuning_param = 1
        
        # Solve complete energy balance 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 8c, Turn on full reaction heat, Time: ", end="")      

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)            
                   
        # ---------------------------------------------------------------------
        # Step 9 - Turn on heat transfer coefficients    
        
        blk.hf.unfix()
        blk.eq_hf.activate()

        blk.hw.unfix()
        blk.eq_hw.activate()

        
        # Solve energy balances w/ film heat transfer 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 9, Turn on gas heat transfer coefficients, Time: ", end="")       

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """)       
            
        # ---------------------------------------------------------------------
        # Step 10 - Turn on film heat transfer between gas and solids              

        # Set value of tuning parameter, to turn down reaction heat
        blk.tuning_param2 = 0.01        
        
        # Unfix gas and solids heat transfer terms and activate constraint  
        blk.Tg_GS.unfix()
        blk.eq_Tg_GS.activate() 
        
        # Solve energy balances w/ film heat transfer 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 10, Turn on 1% gas heat transfer, Time: ", end="")       
            
        results = opt.solve(blk, symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """) 

        # ---------------------------------------------------------------------
        # Step 10a - Turn on heat transfer              

        # Set value of tuning parameter, to turn down reaction heat
        blk.tuning_param2 = 0.1
        
        # Solve energy balances w/ film heat transfer 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 10a, Turn on gas 10% heat transfer, Time: ", end="")         

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """) 
           
        # ---------------------------------------------------------------------
        # Step 10b - Turn on heat transfer              

        # Set value of tuning parameter, to turn down reaction heat
        blk.tuning_param2 = 1
        
        # Solve energy balances w/ film heat transfer 
        ts = time.time()
        if outlvl > 0:
#            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(blk,"- Step 10b, Turn on gas full heat transfer, Time: ", end="")        

        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
        
        if outlvl > 0:
            print(value(time.time() - ts), " s")
            print(results.solver.message)  
            print("""    """) 

        # ---------------------------------------------------------------------
#         Step 11 - Turn on heat loss through refractory-lined reactor wall            
#        
#         Compute overall heat transfer coefficient
#        blk.U.unfix()
#        blk.eq_U.activate()
#        
#        # Heat loss from gas to ambiance 
#        blk.Tg_refractory.unfix()
#        blk.eq_Tg_refractory.activate()        
#        
#        # Solve gas heat energy balance 
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,
#                  "- Step 11, Turn on heat loss to ambiance via refractory...")        
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(results.solver.message)  
#            print("""    """)  
            
# ---------------------------------------------------------------------
        # Step 12 - Turn on wall heat energy balance       

#        # Film heat transfer between gas and reactor wall
#        blk.Tg_GW.unfix()
#        blk.eq_Tg_GW.activate()   
#        
#        # Compute overall heat transfer coefficient (wall-ambiance)
#        blk.Uw.unfix()
#        blk.eq_Uw.activate()
#        
#        # Compute convective gas-wall heat transfer term for wall balance
#        blk.Tw_GW.unfix()
#        blk.eq_Tw_GW.activate()
#
#        # Compute wall-ambiance heat transfer term for wall balance
#        blk.Tw_Wamb.unfix()
#        blk.eq_Tw_Wamb.activate()
#
#        # Solve heat advection 
#        if outlvl > 0:
##            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#            print(blk,"- Step 12, Turn on wall heat balance ...")        
#
#        results = opt.solve(blk,symbolic_solver_labels=False, tee=stee)
#        
#        if outlvl > 0:
#            print(results.solver.message)  
#            print("""    """)  
#
#        # Unfix wall temperature and activate constraint          
#        blk.Tw.unfix()
#        blk.eq_Tw.activate()       
                
