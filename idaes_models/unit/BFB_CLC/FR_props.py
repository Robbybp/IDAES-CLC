# -*- coding: utf-8 -*-
"""
Updated on Tue Dec 20 2018

@authors: cokoli
"""
# Changes the divide behavior to not do integer division
from __future__ import division

# Some more inforation about this module
__author__ = "Chinedu Okoli"
__version__ = "1.0.0"

# Import Pyomo
from pyomo.environ import (Param, Var, Expression, Constraint, 
                Set, Reals, NonNegativeReals, exp, sqrt)

# Import IDAES cores
from idaes_models.core import UnitModel, ProcBlock

@ProcBlock("PropPack")
#class _PropPack(GasProp._PropPack):
class _PropPack(UnitModel):
    """
    This package provides the necessary constraints for solids properties and
    reaction kinetics for oxygen carrier Fe2O3/Fe3O4 on an Al2O3 inert.
    Overall reducer reactions: ## Model assumes all rxns are heterogeneous
    i.e. they all occur in the oxygen carrier
        (1) CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O - Methane combustion
    """
    # Lists of all chemical species and reactions included in package
    scomp = ['Fe2O3','Fe3O4','Al2O3']    
    gcomp = ['CH4','CO2','H2O']
    comp = scomp + gcomp
    rxn_idx = ['1']

    def __init__(self, *args, **kwargs):
        """
        Create a property package object.
        """
        # Extract build arguments
        self.prop_list = kwargs.pop("prop_list",{"gas",
                                                "sol",
                                                "ap",
                                                 "dp",
                                                 "r",
                                                 "h_sol",
                                                 "V_vap",
                                                 "D_vap",
                                                 "k_vap",
                                                 "rho_vap",
                                                 "mu_vap",
                                                 "MW_vap",
                                                 "h_vap"})
        # Call base class constructor
        super(_PropPack,self).__init__(*args,**kwargs)

    def build(self, *args, **kwargs):
        ''' Callable method for Block construction.
        '''
        # Build model
        self._make_params()
        self._make_vars()
        self._make_constraints()
        
    def _make_params(self):
        ''' This section contains a number of parameters or constants
            required by the model.
        '''
        # Set component and reaction liste-2
        self.gcomp = Set(initialize=self.gcomp)
        self.scomp = Set(initialize=self.scomp)
        self.comp = Set(initialize=self.comp)
        self.rxn_idx = Set(initialize=self.rxn_idx)
       
        # Dummy parameter used to scale the rate expression for reaction conversion testing
        self.scale_CH4 = Param(within = NonNegativeReals, default = 1, mutable = True, 
                         doc= 'scale factor for rate rxn of CH4 comb with OC')
        self.eps = Param(default=2e-8,
                      doc='Smoothing Factor for Smooth IF Statements')
        
        # Gas constant
        self.R = Param(default=8.314, doc='Gas constant [J/mol.K]')
        
        self.Tref = Param(default=298.15,
                    doc='Thermodynamic Reference Temperature [K]')
        
        # Mol. weights of solids - units = kg/kmol. ref: NIST webbook
        self.MW = {'CH4':16,'CO':28,'CO2':44,'H2':2,'H2O':18, 
                     'Fe2O3':159.69, 'Fe3O4':231.533, 'Al2O3':101.96}
        
        # Stoichiometric coefficients
        '''Stoichiometric coefficient for each component in each reaction'''
        self.stoic = {('CH4','1'):-1, ('CO2','1'):1,('H2O','1'):2,('Fe2O3','1'):-12,('Fe3O4','1'):8,('Al2O3','1'):0}

        # Std. heat of formation of gas components - units = J/mol - ref: NIST
        self.Hf_comp = {'CH4':-74873.1,'CO2':-393522.4,
                        'H2O':-241826.4,'Fe2O3':-825503.2,'Fe3O4':-1120894,
                        'Al2O3':-1675690}
        
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
    # solid phase only physical xteristics
        self.radg = Param(default = 2.6e-7, mutable = True, 
                          doc = 'rep. particle grain radius within OC part, m') # EPAT
        self.dp = Param(default = 2.8e-4, mutable = True, 
                        doc = 'Diameter of solid particles, m')
        
        # emf - educated guess (EPAT) as rough estimate from ergun equation results (0.4) are suspicious
        self.e_mf = Param(default = 0.45, mutable = True, 
                          doc = 'voidage at minimum fluidization velocity, (-)')
        self.rho_sol = Param(default = 3251.75, mutable = True, 
                             doc = 'Density of solid particles, kg/m^3') #EPAT
        self.k_sol = Param(default = 12.3, mutable = True, 
                             doc = 'Thermal conductivity of solid particles, J/mKs') #EPAT
        
        # vmf - EPAT value used for Davidson model (0.039624)
        self.v_mf = Param(default = 0.039624, mutable = True, 
                          doc = 'velocity at minimum fluidization, m/s')
    # reaction xteristics
        self.rhom = Param(default = 32811, mutable = True, 
                          doc = 'molar density of carrier particle, mol/m^3')
        self.ko = Param(self.rxn_idx, default = {'1':8e-4}, mutable = True, 
                        doc = 'pre-exponential factor, mol^(1-N_rxn)m^(3*N_rxn -2)/s')
        self.E = Param(self.rxn_idx, default = {'1':4.9e4}, mutable = True,
                       doc = 'Activation energy, J/mol')
        self.N_rxn = Param(self.rxn_idx, default = 1.3, mutable = True,
                           doc = 'reaction order in gas species, (-)')
        self.b_rxn = Param(self.rxn_idx, default = 12, mutable = True,
                           doc = 'reaction stoich coeff, (-)')     
        # Available volume for reaction - from EPAT report (1-ep)'
        self.a_vol = Param(within = NonNegativeReals, default = 0.28, mutable = True, 
                         doc= 'available reaction vol. per vol. of OC')

    # gas phase only physical xteristics        
        # Atomic diffusion volumes 
        self.V_diff = {'CH4':24.42,'CO2':26.9,'H2O':13.1} #Ref: (1) Prop gas & liquids (2) Fuller et al. IECR, 58(5), 19, 1966
        
        #Viscosity constants - Reference: Perry and Green Handbook; McGraw Hill, 2008
        self.mu_param = {('CH4','a'):5.2546e-7,('CH4','b'):0.59006,('CH4','c'):105.67,('CH4','d'):0,
                         ('CO2','a'):2.148e-6,('CO2','b'):0.46,('CO2','c'):290,('CO2','d'):0,
                         ('H2O','a'):1.7096e-8,('H2O','b'):1.1146,('H2O','c'):0,('H2O','d'):0}
        
        #Thermal conductivity constants - Reference: Perry and Green Handbook; McGraw Hill, 2008
        self.k_param = {('CH4','a'):8.3983e-6,('CH4','b'):1.4268,('CH4','c'):-49.654,('CH4','d'):0,
                        ('CO2','a'):3.69,('CO2','b'):-0.3838,('CO2','c'):964,('CO2','d'):1.86e6,
                        ('H2O','a'):6.204e-6,('H2O','b'):1.3973,('H2O','c'):0,('H2O','d'):0}

    def _make_vars(self):
        """
        Create property variables.
        """        
        # State variables
        self.P = Var(initialize=101325.0, doc = 'Pressure in Pa') 
        self.y = Var(self.gcomp, initialize=1.0/len(self.gcomp))
        self.x = Var(self.scomp, initialize=1.0/len(self.scomp), bounds = (0,1.6))        
        self.Tg = Var(initialize=298.15, doc = 'Gas Temp in K') 
        self.Ts = Var(initialize=298.15, doc = 'Solid Temp in K')
        # Fixes lower bound of y at 0.5*eps (i.e. 1e-8)
        self.y_max = Var(self.gcomp, initialize=1.0/len(self.gcomp)) 
            
        if "gas" in self.prop_list:
            self.H_comp_g = Var(self.gcomp, domain = Reals, initialize = 1, 
                              doc = 'Gas phase component Enthalpy at Tg(K) - J/mol')            
        if "sol" in self.prop_list:
            self.H_comp_s = Var(self.comp, domain = Reals, initialize = 1, 
                        doc = 'Solid phase component Enthalpy at Ts(K) - J/mol')
                
        # Solid phase only properties
        if "ap" in self.prop_list:
            self.ap = Var(doc = 'Surface area of solid particles, m^2/kg')

        if "h_sol" in self.prop_list:
            self.cp_sol = Var(domain = Reals, doc = 'Heat capacity of solid particles, J/kg.K')
            self.h_sol = Var(within=Reals, initialize=1, doc = 'Enthalpy of solid, J/kg')
            self.Dh_sol = Var(within=Reals, initialize=1, doc = 'Enthalpy change of solid, J/kg')
            self.DH_rxn_s = Var(domain = Reals, initialize = 1, 
                                doc = 'Heat of rxn at system T, J/mol')
            
        # Gas Phase Properties
        if "V_vap" in self.prop_list:
            self.V_vap = Var(domain=NonNegativeReals, initialize = 1,
                             doc = 'gas molar volume - m^3/mol')
        if "D_vap" in self.prop_list:
            self.D_vap = Var(self.gcomp, domain=NonNegativeReals, 
                             doc = 'multicomp diffusion, cm2/s')    
            self.D_bin = Var(self.gcomp, self.gcomp, domain=NonNegativeReals,
                             initialize = 1, doc = 'Binary diffusivity, cm2/s')
            self.temp_7 = Var(self.gcomp, self.gcomp, domain=Reals, 
                               initialize=1, doc = 'dummy var')            
            self.temp_8 = Var(self.gcomp, domain=Reals, 
                               initialize=1, doc = 'dummy var')
            
        if "k_vap" in self.prop_list:
            self.temp_1 = Var(self.gcomp, domain=Reals, 
                               initialize=1e-6, doc = 'dummy var')            
            self.temp_2 = Var(self.gcomp, domain=Reals, 
                               initialize=1e-6, doc = 'dummy var')              
            self.temp_3 = Var(self.gcomp, domain=Reals, 
                               initialize=1e-6, doc = 'dummy var')              
            self.temp_3a = Var(self.gcomp, self.gcomp, domain=NonNegativeReals, 
                               initialize=1, doc = 'dummy var')
            self.k_vap = Var(domain = NonNegativeReals, initialize = 1e-1, 
                             doc = 'thermal conductivity of gas, J/mKs')
            self.k_comp = Var(self.gcomp, domain = NonNegativeReals, initialize = 1e-1, doc = 'pure comp thermal cond, J/mKs')
#            self.A_bin = Var(self.gcomp, self.gcomp, domain=NonNegativeReals,
#                             initialize = 1, doc = 'Binary interaction parameter')
        if "rho_vap" in self.prop_list:
            self.rho_vap = Var(domain=NonNegativeReals, 
                               doc = 'gas mass density - kg/m3')
        if "mu_vap" in self.prop_list:
            self.mu_vap = Var(domain=NonNegativeReals, initialize = 1e1,
                              doc = 'dynamic viscosity of gas, mg/ms')
            self.mu_comp = Var(self.comp, domain=NonNegativeReals, 
                               initialize = 1e1,doc = 'component dynamic viscosity, mg/ms')
            self.temp_4 = Var(self.gcomp, domain=Reals, 
                               initialize=1e-3, doc = 'dummy var')  
            self.temp_5 = Var(self.gcomp, domain=Reals, 
                               initialize=1e-3, doc = 'dummy var')  
            self.temp_6 = Var(self.gcomp, domain=Reals, 
                               initialize=1e-3, doc = 'dummy var') 
        if "MW_vap" in self.prop_list or "rho_vap" in self.prop_list:
            self.MW_vap = Var(domain=Reals, doc = 'mol wt. of gas - kg/mol')
        if "h_vap" in self.prop_list:
            self.cp_vap = Var(domain = NonNegativeReals, 
                              doc = 'heat capacity of gas - J/molK')
            self.h_vap = Var(domain = Reals, doc = 'gas mixture enthalpy - J/mol')
            
        # Reaction Properties
        if "r" in self.prop_list:
            self.X = Var(domain=NonNegativeReals, initialize=0.0, 
                         doc = 'fraction of metal oxide converted')
            self.X_term = Var(domain = Reals, initialize = 1, 
                              doc = 'reformulation term for X to help eqn scaling') 
            self.k = Var(self.rxn_idx, domain=Reals, initialize=1, doc = 'kinetic rate constant')
            self.r_gen = Var(self.rxn_idx, domain=Reals, initialize = 1e-3, 
                             doc = 'gen. rate expression for hetero eqns, units= mol_extent/tonneOC.s')          
            self.rg = Var(self.gcomp, domain=Reals, doc = 'gcomp. total rate expression, units= mol/m^3.s')
            self.rg_rxn = Var(self.gcomp, self.rxn_idx, domain=Reals, doc = 'gcomp. rate expression in each rxn, units= mol/m^3.s')            
            self.rs = Var(self.scomp, domain=Reals, doc = 'scomp. rate expression, units= kg/m^3.s')
            self.C_gas = Var(self.gcomp, domain=NonNegativeReals, initialize=1, doc = 'gas phase molar conc. units = mol/m^3')       
            self.ht_rxn_t = Var(domain = Reals, initialize=1e-6, 
                                doc = 'Total Gas-solids heat transfer due to transfer of reacting mass, J/m^3.s')                        
            self.ht_rxn = Var(self.rxn_idx, domain = Reals, initialize=1e-6, 
                              doc = 'Gas-solids heat transfer due to transfer of reacting mass in a rxn, J/m^3.s')

            self.C_gas_max = Var(self.gcomp, initialize=1) # Fixes lower bound of C_gas at 0.5*eps (i.e. 1e-8)
            self.d1 = Var(self.rxn_idx, domain=Reals, 
                               initialize=1, doc = 'dummy var 1')
            self.d2 = Var(self.rxn_idx, domain=Reals, 
                               initialize=1, doc = 'dummy var 2')
            
    def _make_constraints(self):
        """
        Create property constraints.
        """
        if "gas" in self.prop_list:
            def rule_eq_r7(b,i):
                return self.H_comp_g[i] == 1e3*(b.cp_param[i,'a']*(b.Tg/1000) \
                        + b.cp_param[i,'b']*((b.Tg/1000)**2)/2 \
                        + b.cp_param[i,'c']*((b.Tg/1000)**3)/3 \
                        + b.cp_param[i,'d']*((b.Tg/1000)**4)/4 \
                        - b.cp_param[i,'e']/(b.Tg/1000) + b.cp_param[i,'f'] \
                        - b.cp_param[i,'h'])
            self.eq_r7 = Constraint(self.gcomp, rule=rule_eq_r7)
            
            def rule_eq_q17(b,i):
                return b.cp_param[i,'a'] + b.cp_param[i,'b']*(b.Tg/1000) \
                        + b.cp_param[i,'c']*(b.Tg/1000)**2 \
                        + b.cp_param[i,'d']*(b.Tg/1000)**3 \
                        + b.cp_param[i,'e']/((b.Tg/1000)**2)
            self.cp_comp_g = Expression(self.gcomp, rule=rule_eq_q17)

        # Fixes lower bound of y at 0.5*eps (i.e. 1e-8)      
            def rule_eq_y_lb(b,i):
                return 1e3*b.y_max[i] == 1e3*0.5*(b.y[i] + sqrt(b.y[i]**2 + b.eps**2))
            self.eq_y_lb = Constraint(self.gcomp, rule=rule_eq_y_lb)
            
        if "sol" in self.prop_list:
            def rule_eq_cp_comp_s(b,i):
                return b.cp_param[i,'a'] + b.cp_param[i,'b']*(b.Ts/1000) \
                        + b.cp_param[i,'c']*(b.Ts/1000)**2 \
                        + b.cp_param[i,'d']*(b.Ts/1000)**3 \
                        + b.cp_param[i,'e']/((b.Ts/1000)**2)
            self.cp_comp_s = Expression(self.comp, rule=rule_eq_cp_comp_s) # Solid phase - cp_comp_s units = J/mol.K
        
            def rule_eq_H_comp_s(b,i):
                return self.H_comp_s[i] == 1e3*(b.cp_param[i,'a']*(b.Ts/1000) \
                        + b.cp_param[i,'b']*((b.Ts/1000)**2)/2 \
                        + b.cp_param[i,'c']*((b.Ts/1000)**3)/3 \
                        + b.cp_param[i,'d']*((b.Ts/1000)**4)/4 \
                        - b.cp_param[i,'e']/(b.Ts/1000) + b.cp_param[i,'f'] \
                        - b.cp_param[i,'h'])
            self.eq_H_comp_s = Constraint(self.comp, rule=rule_eq_H_comp_s) # Solid phase - H_comp_s units = J/mol         

        if "ap" in self.prop_list:
            self.eq_ap = Constraint(expr = self.ap*self.rho_sol*self.dp == 6) 

        if "V_vap" in self.prop_list:
            self.eq_q1 = Constraint(expr = self.P*self.V_vap == self.R*self.Tg)
        if "MW_vap" in self.prop_list or "rho_vap" in self.prop_list:
            self.eq_q2 = Constraint(expr = self.MW_vap == 1e-3*sum(self.y_max[i]*self.MW[i] for i in self.gcomp))
        if "rho_vap" in self.prop_list:
            self.eq_q3 = Constraint(expr = self.rho_vap == self.MW_vap*self.P \
                                                        / (self.R*self.Tg))
        if "D_vap" in self.prop_list:          
            def rule_eq_q4(b,i,j):
                return b.D_bin[i,j]*(((b.P*1e-5)*((b.V_diff[i]**(1/3)) \
                        + (b.V_diff[j]**(1/3)))**2)) \
                == (1.43e-3*(b.Tg**1.75) \
                        * ((b.MW[i]+b.MW[j])/(2*b.MW[i]*b.MW[j])) \
                        **0.5) 
            self.eq_q4 = Constraint(self.gcomp, self.gcomp, rule=rule_eq_q4)

            def rule_eq_q5a(b,i,j):
                return 1e4*b.temp_7[i,j]*b.D_bin[i,j] \
                == 1e4*b.y_max[j]
            self.eq_q5a = Constraint(self.gcomp, self.gcomp, rule=rule_eq_q5a)            
            
            def rule_eq_q5(b,i):
                return 1e3*b.D_vap[i]*sum(b.temp_7[i,j] for j in b.gcomp if i != j) \
                == 1e3*(1-b.y_max[i])
            self.eq_q5 = Constraint(self.gcomp, rule=rule_eq_q5)
            
        if "k_vap" in self.prop_list:      
            def rule_eq_q6a(b,i):
                return b.temp_1[i]*b.Tg == b.k_param[i,'c']
            self.eq_q6a = Constraint(self.gcomp, rule=rule_eq_q6a)

            def rule_eq_q6b(b,i):
                return b.temp_2[i]*(b.Tg**2) == b.k_param[i,'d']
            self.eq_q6b = Constraint(self.gcomp, rule=rule_eq_q6b)
            
            def rule_eq_q6c(b,i):
                return b.k_comp[i]*(1 + b.temp_1[i] + b.temp_2[i]) \
                        == b.k_param[i,'a']*(b.Tg**b.k_param[i,'b'])
            self.eq_q6c = Constraint(self.gcomp, rule=rule_eq_q6c) 
            
            def rule_eq_q7b(b,i,j):
                return 1e3*(b.temp_3a[j,i]**2)*b.k_comp[i] == 1e3*b.k_comp[j]
            self.eq_q7b = Constraint(self.gcomp, self.gcomp, rule=rule_eq_q7b)
            
            def rule_eq_q7(b,i,j):
                return (1 + b.temp_3a[j,i] \
                        * ((b.MW[j]/b.MW[i])**0.25))**2 \
                        / (8*(1+(b.MW[j]/b.MW[i])))**0.5    
            self.A_bin = Expression(self.gcomp, self.gcomp, rule=rule_eq_q7)

            def rule_eq_q8b(b,i):
                return 1e3*b.temp_3[i]*(sum(b.y_max[j]*b.A_bin[i,j]**0.5 for j in b.gcomp)) \
                        == 1e3*b.y_max[i]*b.k_comp[i]
            self.eq_q8b = Constraint(self.gcomp, rule=rule_eq_q8b) 
            
            def rule_eq_q8(b):
                return 1e3*self.k_vap == 1e3*sum(b.temp_3[i] for i in b.gcomp)
            self.eq_q8 = Constraint(rule=rule_eq_q8)  

        if "mu_vap" in self.prop_list:            
            def rule_eq_q9a(b,i):
                return b.temp_4[i]*b.Tg == b.mu_param[i,'c']
            self.eq_q9a = Constraint(self.gcomp, rule=rule_eq_q9a)

            def rule_eq_q9b(b,i):
                return b.temp_5[i]*(b.Tg**2) == b.mu_param[i,'d']
            self.eq_q9b = Constraint(self.gcomp, rule=rule_eq_q9b)
            
            def rule_eq_q9c(b,i):
                return b.mu_comp[i]*(1 + b.temp_4[i] + b.temp_5[i]) \
                        == (1e6*b.mu_param[i,'a'])*(b.Tg**b.mu_param[i,'b'])
            self.eq_q9c = Constraint(self.gcomp, rule=rule_eq_q9c)             

            def rule_eq_q10b(b,i):
                return b.temp_6[i]*(sum(b.y_max[j]*(b.MW[j]/b.MW[i])**0.5 for j in b.gcomp)) \
                        == b.y_max[i]*b.mu_comp[i]
            self.eq_q10b = Constraint(self.gcomp, rule=rule_eq_q10b) 

            def rule_eq_q10(b):
                return self.mu_vap == sum(b.temp_6[i] for i in b.gcomp)
            self.eq_q10 = Constraint(rule=rule_eq_q10)

        if "r" in self.prop_list:             
            def rule_eq_q13(b, i):
                return b.C_gas[i]*b.R*b.Tg == b.y[i]*b.P
            self.eq_q13 = Constraint(self.gcomp, rule=rule_eq_q13, doc = 'gas phase conc. calc.') 

            # Fixes lower bound of C_gas at 0.5*eps (i.e. 1e-8)        
            def rule_eq_C_gas_lb(b,i):
                return b.C_gas_max[i]*1e3 == 1e3*0.5*(b.C_gas[i] + sqrt(b.C_gas[i]**2 + b.eps**2))
            self.eq_C_gas_lb = Constraint(self.gcomp, rule=rule_eq_C_gas_lb)
                                                                                          
            # rate expression constraints - overall rate and species specific rates
            def rule_eq_r1(b, i):
                return 1e6*b.k[i] == 1e6*b.ko[i]*b.d1[i]
            self.eq_r1 = Constraint(self.rxn_idx, rule=rule_eq_r1, 
                                    doc = 'kinetic rate constant eqn')  
            def rule_eq_r1b(b, i):
                return 1e3*b.d1[i] == 1e3*exp(b.d2[i])
            self.eq_r1b = Constraint(self.rxn_idx, rule=rule_eq_r1b, doc = '')             
            def rule_eq_r1c(b,i):
                return 1e-3*b.d2[i]*b.R*b.Ts == 1e-3*-b.E[i]
            self.eq_r1c = Constraint(self.rxn_idx, rule=rule_eq_r1c, doc = '')           


            """
            This equation is for calculating the conv. of fe2 to fe3 in the reactor. It is gotten
            by converting the EPAT conv. equation (based on mole flow) to a mass fraction equivalent.
            
            X =                         x_fe3
                ---------------------------------------------------------
                 x_fe3 + x_fe2*((stoich_fe3/-stoich_fe2)*(MW_fe3/MW_fe2))
                
            where  x_fe2,x_fe3 = wt. fraction of fe2 and fe3 at any
                    given point in the reactor (state variables)                     
            """  

            def rule_eq_r2(b):
                return 1e3*b.X*(b.x['Fe3O4'] + (b.MW['Fe3O4']/b.MW['Fe2O3']) \
                                *(b.stoic['Fe3O4','1']/-b.stoic['Fe2O3','1'])*b.x['Fe2O3']) \
                        == 1e3*b.x['Fe3O4']
            self.eq_r2 = Constraint(rule=rule_eq_r2, doc = 'conversion of metal oxide eqn')  

            def rule_eq_r2b(b):
                return 1e6*b.X_term**3 == 1e6*(1-b.X)**2   
            self.eq_r2b = Constraint(rule = rule_eq_r2b)

           
            # general reaction rate expression of CH4 with Fe2O3
            def rule_eq_r4(b,i):
                return 1e3*b.r_gen[i]== 1e3*1e6*b.scale_CH4*b.x['Fe2O3']*(b.a_vol/b.MW['Fe2O3'])\
                                    *3*b.b_rxn[i]*b.k[i]*((b.C_gas_max['CH4']+1e-8)**b.N_rxn[i])\
                                    *b.X_term/(b.rhom*b.radg)/-b.stoic['Fe2O3','1']                   
            self.eq_r4 = Constraint(self.rxn_idx, rule=rule_eq_r4, 
                                    doc = 'general rate expression (mole extent basis), \
                                            units are mol_extent/tonneOC.s')        

            def rule_eq_r5(b,i):
                return 1e3*b.rs[i] == 1e3*b.rho_sol*b.MW[i] \
                                    *sum(b.stoic[i,j]*1e-6*b.r_gen[j] for j in b.rxn_idx)
            self.eq_r5 = Constraint(self.scomp, rule=rule_eq_r5, 
                                    doc = 'comp specific rate expression')

            def rule_eq_r6(b,i,j):
                return 1e3*b.rg_rxn[i,j] == 1e3*-b.stoic[i,j]*(1e3*b.rho_sol*1e-6*b.r_gen[j])                                 
            self.eq_r6 = Constraint(self.gcomp, self.rxn_idx, rule=rule_eq_r6, 
                                    doc = 'comp specific rate expression for a rxn')

            def rule_eq_r6b(b,i):
                return 1e3*b.rg[i] == 1e3*sum(-b.stoic[i,j]*(1e3*b.rho_sol*1e-6*b.r_gen[j]) \
                                               for j in b.rxn_idx)
            self.eq_r6b = Constraint(self.gcomp, rule=rule_eq_r6b,
                                     doc = 'total comp specific rate expression') 

            def rule_eq_r8(b,i):
                return b.ht_rxn[i] == b.rg_rxn['CH4',i]*b.H_comp_g['CH4'] \
                        + b.rg_rxn['H2O',i]*b.H_comp_s['H2O'] \
                        + b.rg_rxn['CO2',i]*b.H_comp_s['CO2'] 
            self.eq_r8 = Constraint(self.rxn_idx, rule = rule_eq_r8)

            self.eq_r8a = Constraint(expr = self.ht_rxn_t \
                                     == sum(self.ht_rxn[i] for i in self.rxn_idx))            
            
        if "h_vap" in self.prop_list:
#            self.cp_vap.fix(59.2189)   
            self.eq_q18 = Constraint(expr = self.cp_vap == \
                                    sum(self.cp_comp_g[j]*self.y_max[j] for j in self.gcomp))             
             

            self.eq_q21 = Constraint(expr = self.h_vap \
                                        == sum((self.H_comp_g[j])*self.y[j] for j in self.gcomp))            
            
        if "h_sol" in self.prop_list:                    
            self.eq_h1 = Constraint(expr = self.cp_sol \
                                    == sum(self.cp_comp_s[j]*(1000/self.MW[j]) \
                                           *self.x[j] for j in self.scomp))   
            
            self.eq_h2 = Constraint(expr = 1e-3*self.Dh_sol \
                                    == 1e-3*sum(self.H_comp_s[j]*(1000/self.MW[j]) \
                                           *self.x[j] for j in self.scomp))             

            def rule_eq_h3(b):
                return b.DH_rxn_s == sum(b.stoic[i,'1']*(b.H_comp_s[i] + b.Hf_comp[i]) for i in b.comp)
            self.eq_h3 = Constraint(rule = rule_eq_h3)  
                                      
            self.eq_h4 = Constraint(expr = self.h_sol == self.Dh_sol + self.DH_rxn_s*1000 \
                                    *(self.x['Fe3O4']/(self.stoic['Fe3O4','1']*self.MW['Fe3O4']))) 