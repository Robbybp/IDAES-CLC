#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:55:31 2019

A simple flowsheet model for the simulation of a methane-fueled MB fuel 
reactor.

@author: Anca Ostace (aostace)
"""

from __future__ import division
from __future__ import print_function

__author__ = "Chinedu Okoli and Anca Ostace"
__version__ = "2.0.0"

from pyomo.environ import value, Var, Constraint
from pyomo.environ import ConcreteModel
from pyomo.environ import Block
from pyomo.core.base.sets import _SetProduct
from pyomo.core.base.constraint import SimpleConstraint
from pyomo.core.base.var import SimpleVar
from pyomo.core.base.indexed_component import IndexedComponent
#from pyomo.environ import *
from pyomo.opt import SolverFactory

import time

from idaes_models.core import FlowsheetModel, ProcBlock
import MB_CLC as MB_CLC_fuel
import ss_sim
from CLC_integrate import alg_update, integrate

import pdb

@ProcBlock("Flowsheet")
class _Flowsheet(FlowsheetModel):
    def __init__(self, *args, **kwargs):
        """
        Create a flowsheet model.
        """
        FlowsheetModel.__init__(self, *args, **kwargs)

    def build(self):
        """
        Make the flowsheet object, fix some variables, and solve the problem
        """

        # Create a custom grid, fe_set        
        # why create an fe_set instead of using Transformation 
        # factory?
        nfe = 6
        fe_a = 1/4.0
        fe_b = 0.2
        fe_set = [0, 0.004]
        for i in range(1,nfe+1):
            if i < nfe*fe_a:
                fe_set.append(i*fe_b/(nfe*fe_a))
            elif i == nfe:   
                fe_set.append(1)
            else:
                fe_set.append(fe_b + (i-nfe*fe_a)*(1-fe_b)/(nfe*(1-fe_a)))

        """
        Args: (to MB_CLC_fuel object, as defined in model file)
        dae_method = method to use for calcuating derivatives (default = OCLR)
                    - BFD1 - 1st order backwards finite difference
                    - OCLR - Orthogonal collocation, Lagrange-Radau
                    - OCLL - Orthogonal collocation, Lagrange-Legendre
        press_drop = Pressure drop correlation for superficial velocity calc.
                    - SimplifiedP - simplified pressure correlations 
                    - Ergun - Ergun equation
        fe_set = set of normalised finite element locations
        nfe = number of finite elements for bed discretization (default = 15)
                    (not used if fe_set specified)
        ncp = number of collocation points (OCLR or OCLL only, default = 3)

        fe_set_t
        nfe_t 
        ^ adding time set as a model-level continuous set...
        will change when moving to new framework
        """           
        # Create unit model for fuel reactor
        # unit model - an attribute of the flowsheet model
        # "a block within a block"
        #self.MB_fuel = MB_CLC_fuel.MB(
        #        parent=self,
        #        dae_method = 'OCLR',
        #        press_drop = 'Ergun',
        #        fe_set = fe_set,
        #        ncp = 3)

        # need to add time set to the above
        # open question still:
        # how long of a horizon should I simulate?
        #
        # why is nfe in z-dimension not an input here?
        self.MB_fuel = MB_CLC_fuel.MB(
                parent=self,
                dae_method = 'OCLR',
                press_drop = 'Ergun',
                fe_set = fe_set,
                ncp = 3,
                horizon = 1, # was 10
                nfe_t = 1,   #  "  "
                ncp_t = 3)


        
def setInputs(fs):    
    # ===== Fuel Reactor ===== 
    # Gas phase inlet conditions
    for t in fs.MB_fuel.t:
        fs.MB_fuel.Gas_In_F[t].fix(128.20513)   # mol/s
        fs.MB_fuel.Gas_In_P[t].fix(2.00)   # bar *estimated. min pressure to overcome pressure drop
        fs.MB_fuel.Gas_In_Tg[t].fix(293.15)     # K
        fs.MB_fuel.Gas_In_y['CO2',t].fix(0.02499)
        fs.MB_fuel.Gas_In_y['H2O',t].fix(0.00001)
        fs.MB_fuel.Gas_In_y['CH4',t].fix(0.975) 
           
        # Solid phase inlet conditions
        fs.MB_fuel.Solid_In_M[t].fix(591.4) #479.011) # kg/s
        fs.MB_fuel.Solid_In_Ts[t].fix(1183.15)      # K
        fs.MB_fuel.Solid_In_x['Fe2O3',t].fix(0.45)
        fs.MB_fuel.Solid_In_x['Fe3O4',t].fix(1e-9)
        fs.MB_fuel.Solid_In_x['Al2O3',t].fix(0.55)
        
    # Bed characteristics
    fs.MB_fuel.Dr.fix(6.5) # m
    fs.MB_fuel.L.fix(5) # m
    fs.MB_fuel.eps.fix(0.4) # (-)   

def perturbInputs(fs):
    t_1 = fs.MB_fuel.t.get_finite_elements()[1]
    for t in fs.MB_fuel.t:
        if t < t_1:
            #print(t)
            fs.MB_fuel.Solid_In_M[t].fix(691.4)
            #fs.MB_fuel.Solid_In_x['Fe2O3',t].fix(0.46)
            #fs.MB_fuel.Solid_In_x['Al2O3',t].fix(0.54)

def setICs(fs,fs_ss):
    # getting the names from the variables would only be useful if I have a set
    # of differential variables defined already
    diff_vars_t = []
    diff_vars_t.append('Cg')
    diff_vars_t.append('q')
    diff_vars_t.append('Tg')
    diff_vars_t.append('Ts')

    for var_ss in fs_ss.MB_fuel.component_objects(Var,active=True):
        var_name = var_ss.getname()
        if var_name in diff_vars_t:

            if type(var_ss.index_set()) is _SetProduct:
                ss_index_sets = var_ss.index_set().set_tuple
            else: 
                ss_index_sets = var_ss.index_set()

            ic_param = getattr(fs.MB_fuel,var_name+'_0')

            print(ic_param)
            for index in var_ss:

                if index is None:
                    ss_value = value(var_ss)
                    ic_param.set_value(ss_value)
                else:
                    ss_value = value(var_ss[index])
                ic_param[index].set_value(ss_value)

def initialize_ss(fs,fs_ss):

    time = fs.MB_fuel.t
    for var_ss in fs_ss.MB_fuel.component_objects(Var,active=True):

        var_name = var_ss.getname()
        var_ol = getattr(fs.MB_fuel,var_name)

        if type(var_ol.index_set()) is _SetProduct:
            ol_index_sets = var_ol.index_set().set_tuple
        else:
            ol_index_sets = var_ol.index_set()
        # ^ can end up being a tuple of sets or just a pyomo set
        # if ol var is not indexed, this guy is {None}

        # load value of the ss variable, for each ss index, into the 
        # appropriate open loop variable
        #pdb.set_trace()
        for index in var_ss:
        # for an unindexed variable, this is [None]
        # and it seems the loop is skipped...
        # (maybe this is for fixed variables)
        # ^ better way to do this: check if var is SimpleVar
            if var_ss[index].stale == False:
                ss_value = value(var_ss[index])
            else:
                continue
            index_type = type(index)

            if index is None:
                if time in ol_index_sets:
                    for t in time:
                        var_ol[t].set_value(ss_value)
                else: 
                    var_ol.set_value(ss_value)
                continue

            elif index_type is tuple:
                if time in ol_index_sets:
                    for t in time: 
                        ol_index = index + (t,)
                        var_ol[ol_index].set_value(ss_value)
                else:
                    var_ol[index].set_value(ss_value)
                continue

            # here, really want to check if ss_var is indexed by a single set
            # doesn't matter what type that is 
            # so should I check if index_type is not tuple? probably
            # (just 'else' would be fine)

            #elif index_type is int or index_type is float:
            else:
                if time in ol_index_sets:
                    for t in time:
                        ol_index = (index,t)
                        var_ol[ol_index].set_value(ss_value)
                else: 
                    var_ol[index].set_value(ss_value)
                continue

#def alg_update(fs,t):



        
def print_summary_fuel_reactor(fs):
    """
    Print some key results.  """
    print("\nResults:")
    print("==========================================")
    print("---Moving Bed Fuel Reactor---")       
    
    print("\nInlet gas: ", 
              "\nCO2: ", value(fs.MB_fuel.F[0,'CO2',0]), "mol/s",
              "\nH20: ", value(fs.MB_fuel.F[0,'H2O',0]), "mol/s",
              "\nCH4: ", value(fs.MB_fuel.F[0,'CH4',0]), "mol/s",
              "\nCO2: ", value(fs.MB_fuel.Gas_M[0,'CO2',0]), "kg/s",
              "\nH20: ", value(fs.MB_fuel.Gas_M[0,'H2O',0]), "kg/s",
              "\nCH4: ", value(fs.MB_fuel.Gas_M[0,'CH4',0]), "kg/s")
    print("\nOutlet gas: ", 
              "\nCO2: ", value(fs.MB_fuel.F[1,'CO2',0]), "mol/s",
              "\nH20: ", value(fs.MB_fuel.F[1,'H2O',0]), "mol/s", 
              "\nCH4: ", value(fs.MB_fuel.F[1,'CH4',0]), "mol/s",
              "\nCO2: ", value(fs.MB_fuel.Gas_M[1,'CO2',0]), "kg/s",
              "\nH20: ", value(fs.MB_fuel.Gas_M[1,'H2O',0]), "kg/s", 
              "\nCH4: ", value(fs.MB_fuel.Gas_M[1,'CH4',0]), "kg/s")
    print("\nInlet solids: ", 
              "\nFe2O3: ", value(fs.MB_fuel.Solid_F[1,'Fe2O3',0]), "mol/s",
              "\nFe3O4: ", value(fs.MB_fuel.Solid_F[1,'Fe3O4',0]), "mol/s", 
              "\nAl: ", value(fs.MB_fuel.Solid_F[1,'Al2O3',0]), "mol/s",
              "\nFe2O3: ", value(fs.MB_fuel.Solid_M[1,'Fe2O3',0]), "kg/s",
              "\nFe3O4: ", value(fs.MB_fuel.Solid_M[1,'Fe3O4',0]), "kg/s", 
              "\nAl: ", value(fs.MB_fuel.Solid_M[1,'Al2O3',0]), "kg/s")
    print("\nOutlet solids: ", 
              "\nFe2O3: ", value(fs.MB_fuel.Solid_F[0,'Fe2O3',0]), "mol/s",
              "\nFe3O4: ", value(fs.MB_fuel.Solid_F[0,'Fe3O4',0]), "mol/s", 
              "\nAl: ", value(fs.MB_fuel.Solid_F[0,'Al2O3',0]), "mol/s",
              "\nFe2O3: ", value(fs.MB_fuel.Solid_M[0,'Fe2O3',0]), "kg/s",
              "\nFe3O4: ", value(fs.MB_fuel.Solid_M[0,'Fe3O4',0]), "kg/s", 
              "\nAl: ", value(fs.MB_fuel.Solid_M[0,'Al2O3',0]), "kg/s") 
    
    print("\nGas inlet velocity: ", value(fs.MB_fuel.vg[0,0]), "m/s")
    print("Gas outlet velocity: ", value(fs.MB_fuel.vg[1,0]), "m/s")
    print("Solids velocity: ", value(fs.MB_fuel.vs[0]), "m/s")    
    
    print("\nHeat of reaction @ z=0: ", 
              value(fs.MB_fuel.DH_rxn_s[0,0]), "J/(mol reaction)")
    print("Heat of reaction @ z=1: ", 
              value(fs.MB_fuel.DH_rxn_s[1,0]), "J/(mol reaction)")
    
    print("\nCH4 conversion: ", value(fs.MB_fuel.X_gas[0])*100, " %")
    print("Fe2O3 conversion: ", value(fs.MB_fuel.X_OC[0])*100, " %")
    
    print('\nPressure @inlet: ', value(fs.MB_fuel.P[0,0]))
    print('Pressure @outlet: ', value(fs.MB_fuel.Gas_Out_P[0]))
    
    print("\nReactor bed height:", value(fs.MB_fuel.L), " m")
    print("Reactor bed diameter:", value(fs.MB_fuel.Dr), " m")
#    print("Refractory wall thickness", value(fs.MB.refractory_th), " m")
    
    print("\nInlet gas flow:", value(fs.MB_fuel.Gas_In_F[0]), " mol/s")
    print("Outlet gas flow:", value(fs.MB_fuel.Ftotal[1,0]), " mol/s")
    print("Inlet solids flow:", value(fs.MB_fuel.Solid_In_M[0]), " kg/s")
    print("Outlet solids flow:", value(fs.MB_fuel.Solid_Out_M[0]), " kg/s")
    print("Inlet solids temperature:", value(fs.MB_fuel.Solid_In_Ts[0]), " K")
    print("Outlet solids temperature:", value(fs.MB_fuel.Solid_Out_Ts[0]), " K")
    
    print("Inlet gas temperature:", value(fs.MB_fuel.Tg[0,0]), " K")
    print("Outlet gas temperature:", value(fs.MB_fuel.Tg[1,0]), " K")    
    
    print("\nInlet solid mass fractions: ", 
              "\nFe2O3: ", value(fs.MB_fuel.x[1,'Fe2O3',0]),
              "\nFe3O4: ", value(fs.MB_fuel.x[1,'Fe3O4',0]), 
              "\nAl2O3: ", value(fs.MB_fuel.x[1,'Al2O3',0]))
    print("Outlet solid mass fractions: ", 
              "\nFe2O3: ", value(fs.MB_fuel.x[0,'Fe2O3',0]),
              "\nFe3O4: ", value(fs.MB_fuel.x[0,'Fe3O4',0]), 
              "\nAl2O3: ", value(fs.MB_fuel.x[0,'Al2O3',0])) 
    
#
#def results_plot_fuel_reactor(self):
#    """
#    Plot some key results.
#    """
#    
#    import matplotlib.pyplot as plt    
#
#    # Total pressure profile
#    P = []
#    for z in self.MB_fuel.z:
#        P.append(value(self.MB_fuel.P[z]))
#    fig_P = plt.figure(1)
#    plt.plot(self.MB_fuel.z, P)
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Total Pressure [bar]")       
#
#    # Temperature profile
#    Tg = []
#    Ts = []
##    Tw = []
#    for z in self.MB_fuel.z:
#        Tg.append(value(self.MB_fuel.Tg[z] - 273.15))
#        Ts.append(value(self.MB_fuel.Ts[z] - 273.15))
##        Tw.append(value(self.MB_fuel.Tw[z]))
#    fig_T = plt.figure(2)
#    plt.plot(self.MB_fuel.z, Tg, label='Tg')
#    plt.plot(self.MB_fuel.z, Ts, label='Ts')
##    plt.plot(self.MB_fuel.z, Tw, label='Tw')
#    plt.legend(loc=0,ncol=2)
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Temperature [C]") 
#    
#    # Superficial gas velocity and minimum fluidization velocity
#    vg = []
#    umf = []
#    for z in self.MB_fuel.z:
#        vg.append(value(self.MB_fuel.vg[z]))
#        umf.append(value(self.MB_fuel.umf[z]))
#    fig_vg = plt.figure(3)
#    plt.plot(self.MB_fuel.z, vg, label='vg')
#    plt.plot(self.MB_fuel.z, umf, label='umf')
#    plt.legend(loc=0,ncol=2)
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Superficial gas velocity [m/s]")
#    
#   # Gas components molar flow rate
#    for j in self.MB_fuel.GasList:
#        F = []
#        for z in self.MB_fuel.z:
#            F.append(value(self.MB_fuel.F[z,j]))
#        fig_F = plt.figure(4)
#        plt.plot(self.MB_fuel.z, F, label=j)
#    plt.legend(loc=0,ncol=len(self.MB_fuel.GasList))
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Gas component molar flow rate, F [mol/s]")  
#    
#    # Bulk gas phase total molar flow rate
#    Ftotal = []
#    for z in self.MB_fuel.z:
#        Ftotal.append(value(self.MB_fuel.Ftotal[z]))
#    fig_Ftotal = plt.figure(5)
#    plt.plot(self.MB_fuel.z, Ftotal)
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Total molar gas flow rate [mol/s]")  
#
#    # Solid components mass flow rate
#    for j in self.MB_fuel.SolidList:
#        M = []
#        for z in self.MB_fuel.z:
#            M.append(value(self.MB_fuel.Solid_M[z,j]))
#        fig_M = plt.figure(6)
#        plt.plot(self.MB_fuel.z, M, label=j)
#    plt.legend(loc=0,ncol=len(self.MB_fuel.SolidList))
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Solid components mass flow rate [kg/s]")
#    
#     # Bulk solid phase total molar flow rate
#    Mtotal = []
#    for z in self.MB_fuel.z:
#        Mtotal.append(value(self.MB_fuel.Solid_M_total[z]))
#    fig_Mtotal = plt.figure(7)
#    plt.plot(self.MB_fuel.z, Mtotal)
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Solid total mass flow rate [kg/s]")        
#    
#    # Gas phase concentrations
#    for j in self.MB_fuel.GasList:
#        Cg = []
#        for z in self.MB_fuel.z:
#            Cg.append(value(self.MB_fuel.Cg[z,j]))
#        fig_Cg = plt.figure(8)
#        plt.plot(self.MB_fuel.z, Cg, label=j)
#    plt.legend(loc=0,ncol=len(self.MB_fuel.GasList))
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Concentration [mol/m3]")       
#    
#    # Gas phase mole fractions
#    for j in self.MB_fuel.GasList:
#        y = []
#        for z in self.MB_fuel.z:
#            y.append(value(self.MB_fuel.y[z,j]))
#        fig_y = plt.figure(9)
#        plt.plot(self.MB_fuel.z, y, label=j)
#    plt.legend(loc=0,ncol=len(self.MB_fuel.GasList))
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("y [-]")  
#    
#    # Solid phase mass fractions
#    for j in self.MB_fuel.SolidList:
#        x = []
#        for z in self.MB_fuel.z:
#            x.append(value(self.MB_fuel.x[z,j]))
#        fig_x = plt.figure(10)
#        plt.plot(self.MB_fuel.z, x, label=j)
#    plt.legend(loc=0,ncol=len(self.MB_fuel.SolidList))
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("x [-]")  
#
#    # Total mass fraction
#    xtot = []
#    for z in self.MB_fuel.z:
#        xtot.append(value(self.MB_fuel.xtot[z]))
#    fig_xtot = plt.figure(11)
#    plt.plot(self.MB_fuel.z, xtot)
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Total mass fraction [-]") 
#    
#    # # Gas mix density
#    # rhog = []
#    # for z in self.MB_fuel.z:
#    #     rhog.append(value(self.MB_fuel.rho_vap[z]))
#    # fig_rhog = plt.figure(23)
#    # plt.plot(self.MB_fuel.z, rhog)
#    # plt.grid()
#    # plt.xlabel("Bed height [-]")
#    # plt.ylabel("Gas mix density [kg/m3]") 
#               
#    # Fe conversion
#    X_Fe = []
#    for z in self.MB_fuel.z:
#        X_Fe.append(value(self.MB_fuel.X[z])*100)
#    fig_X_Fe = plt.figure(13)
#    plt.plot(self.MB_fuel.z, X_Fe)
#    plt.grid()
#    plt.xlabel("Bed height [-]")
#    plt.ylabel("Fraction of metal oxide converted [%]") 
#                       
       
def main():
    """
    Make the flowsheet object and solve
    """
    ss_flowsheet = ss_sim.main()
    flowsheet = Flowsheet(name='MB_Model') 

    # fill in values of IC parameters from steady state solve
    setICs(flowsheet,ss_flowsheet)
    
    # Fix variables
    setInputs(flowsheet) 

    # Initialize at steady state
    initialize_ss(flowsheet,ss_flowsheet)

    mb = flowsheet.MB_fuel

    # Then perturb
    # ^ that function should go in this file, probably
    # input: dict mapping state names (strings) to new values

    # this seems like as much work as doing the perturbation...
    # maybe just make a self contained function
    #input_perturbation = { ('Solid_In_x',mb.t) : { ['Fe2O3'] : 0.25 },
    #                       ('Solid_In_x',mb.t) : { ['Al2O3'] : 0.75 } }

    perturbInputs(flowsheet)

    # perturb states

    # should put this in a dedicated ~intialize~ function
    # that also intelligently initializes the model after perturbation
    mb.eq_d4.deactivate()
    mb.eq_d5.deactivate()
    mb.eq_d8.deactivate()
    mb.eq_d9.deactivate()
    mb.eq_d10.deactivate()
    mb.eq_g7.deactivate()
    mb.eq_g8.deactivate()
    mb.eq_g10.deactivate()
    mb.eq_g11.deactivate()
    mb.eq_g12.deactivate()
    mb.eq_g13.deactivate()
    mb.eq_g14.deactivate()
    mb.eq_g4.deactivate()
    mb.eq_g5.deactivate()
    mb.eq_g2.deactivate()
    mb.Tg_GW.fix(0.0)
    mb.Tg_refractory.fix(0.0)
    mb.Tw_Wamb.fix()
    mb.Tw.fix()
    mb.Nuw.fix()
    mb.Nu_ext.fix()
    mb.hw.fix()
    mb.hext.fix()
    mb.hext2.fix()
    mb.U.fix()
    mb.Uw.fix()
    mb.Pr_ext.fix()
    mb.Ra.fix()
    mb.Re.fix()
    ###
    

    '''
    ts = time.time() 
    # Initialize fuel reactor
    flowsheet.MB_fuel._initialize(outlvl=1,
                              optarg={"tol"            : 1e-8,
                                      "max_cpu_time"   : 600,
                                      "print_level"    : 5,
                                      "halt_on_ampl_error": 'yes'})        
    '''

    # Create a solver
    opt = SolverFactory('ipopt')
    opt.options = {'tol': 1e-8,
                   'linear_solver' : 'ma27',
                   'bound_push': 1e-8,
                   'max_cpu_time': 600,
                   'print_level': 5,
                   'halt_on_ampl_error': 'yes'}
    flowsheet.write('fs.nl')

    with open('dyn_fs.txt','w') as f:
        flowsheet.display(ostream=f)

    print('Constraints violated pre-solve:')
    for const in flowsheet.MB_fuel.component_objects(Constraint,active=True):
        if not isinstance(const,SimpleConstraint):
            for idx in const:
                if (value(const[idx].body) > value(const[idx].upper) + 1.0e-7) or \
                        (value(const[idx].body) < value(const[idx].lower) - 1.0e-7):
                    print(const.name,idx)
        else:
            if (value(const.body) > value(const.upper) + 1.0e-7) or \
                    (value(const.body) < value(const.lower) - 1.0e-7):
                print(const.name)
    print('- - -\n')

    print('Variable bounds violated pre-solve:')
    for var in flowsheet.MB_fuel.component_objects(Var,active=True):
        # don't use IndexedComponent here, variables are always indexed components 
        # could also do this by iterating over component_objects(SimpleVar)...?
        if not isinstance(var,SimpleVar):
            for idx in var:
                if not (var[idx].lb is None):
                    if (var[idx].value < var[idx].lb - 1.0e-7):
                        pdb.set_trace()
                        print(var.name,idx)
                if not (var[idx].ub is None):
                    if (var[idx].value > var[idx].ub + 1.0e-7):
                        pdb.set_trace()
                        print(var.name,idx)
        else:
            if (var.value > var.ub + 1.0e-7) or \
                    (var.value < var.lb - 1.0e-7):
                print(var.name)
    print('- - -\n')

    # initialized at steady state, works regardless:
    flowsheet.strip_bounds()

    for z in mb.z:
        for t in mb.t:
            mb.Cg[z,'CH4',t].setlb(1e-6)


    # want a function to integrate the model one step from a specified time point
    # will call for all t
    integrate(flowsheet,0)
    
    #results = opt.solve(flowsheet,tee=True,symbolic_solver_labels=False,
    #                        keepfiles=False)

    #with open('dyn_fs_sol.txt','w') as f:
    #    flowsheet.display(ostream=f)


    
    '''
    
    print("\n")
    print("----------------------------------------------------------")
    print('Total simulation time: ', value(time.time() - ts), " s")
    print("----------------------------------------------------------")

    
    # Print some variables 
    print_summary_fuel_reactor(flowsheet) 

    # Plot some variables 
    #results_plot_fuel_reactor(flowsheet) 

    #with open('m_fs.txt','w') as f:
    #    flowsheet.display(ostream=f)

    # Store the flowsheet    
    '''
    return flowsheet
    
if __name__ == "__main__":
    flowsheet = main()         
