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

from pyomo.environ import value, Var, Constraint, Objective
from pyomo.environ import ConcreteModel
from pyomo.environ import Block
from pyomo.core.base.sets import _SetProduct
from pyomo.core.base.constraint import SimpleConstraint
from pyomo.core.base.var import SimpleVar
from pyomo.core.base.indexed_component import IndexedComponent
#from pyomo.environ import *
from pyomo.opt import SolverFactory

from pyomo.contrib.pynumero.sparse import BlockSymMatrix
from pyomo.contrib.pynumero.interfaces.nlp import PyomoNLP
from pyomo.contrib.pynumero.extensions import hsl

from pyomo.dae import DerivativeVar, ContinuousSet, Simulator

import time

import scipy as sp
import numpy as np
import casadi
import matplotlib.pyplot as plt

from idaes_models.core import FlowsheetModel, ProcBlock
import MB_CLC as MB_CLC_fuel
import ss_sim
from CLC_integrate import alg_update, integrate, update_time_derivatives, implicit_integrate

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
        # controlled by fe_set...
        self.MB_fuel = MB_CLC_fuel.MB(
                parent=self,
                dae_method = 'OCLR',
                press_drop = 'Ergun',
                fe_set = fe_set,
                ncp = 3,
                horizon = 10, # was 10, then 1, then 10^-2, then 10^-4, now back to 1...
                nfe_t = 10,   #  "  "
                ncp_t = 3) # was 3

        
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
        fs.MB_fuel.Solid_In_x['Fe2O3',t].fix(0.44999)
        fs.MB_fuel.Solid_In_x['Fe3O4',t].fix(1e-5)
        fs.MB_fuel.Solid_In_x['Al2O3',t].fix(0.55)
        
    # Bed characteristics
    fs.MB_fuel.Dr.fix(6.5) # m
    fs.MB_fuel.L.fix(5) # m
    fs.MB_fuel.eps.fix(0.4) # (-)   

def perturbInputs(fs,t,**kwargs):
    m = fs.MB_fuel
    if 'Solid_M' in kwargs:
        m.Solid_In_M[t].fix( kwargs['Solid_M'] )
    if 'Solid_T' in kwargs:
        m.Solid_In_Ts[t].fix( kwargs['Solid_T'] )
    if 'Solid_x' in kwargs:
        m.Solid_In_x['Fe2O3',t].fix( kwargs['Solid_x']['Fe2O3'] )
        m.Solid_In_x['Fe3O4',t].fix( kwargs['Solid_x']['Fe3O4'] )
        m.Solid_In_x['Al2O3',t].fix( kwargs['Solid_x']['Al2O3'] )
    if 'Gas_F' in kwargs:
        m.Gas_In_F[t].fix( kwargs['Gas_F'] ) 
    if 'Gas_P' in kwargs:
        m.Gas_In_P[t].fix( kwargs['Gas_P'] )
    if 'Gas_T' in kwargs:
        m.Gas_In_Tg[t].fix( kwargs['Gas_T'] )
    if 'Gas_y' in kwargs:
        m.Gas_In_y['CO2',t].fix( kwargs['Gas_y']['CO2'] )
        m.Gas_In_y['H2O',t].fix( kwargs['Gas_y']['H2O'] )
        m.Gas_In_y['CH4',t].fix( kwargs['Gas_y']['CH4'] )

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
    
def print_violated_constraints(flowsheet,tol=1.0e-8):

    print('\nConstraints violated:')
    for const in flowsheet.MB_fuel.component_objects(Constraint,active=True):
        if not isinstance(const,SimpleConstraint):
            for idx in const:
                up_infeas = value(const[idx].upper) - value(const[idx].body)
                lo_infeas = value(const[idx].body) - value(const[idx].lower)
                if (value(const[idx].body) > value(const[idx].upper) + tol) or \
                        (value(const[idx].body) < value(const[idx].lower) - tol):
                    print(const.name,idx,value(const[idx].body))
        else:
            if (value(const.body) > value(const.upper) + tol) or \
                    (value(const.body) < value(const.lower) - tol):
                print(const.name)
    print('- - -\n')

    print('Variable bounds violated:')
    for var in flowsheet.MB_fuel.component_objects(Var,active=True):
        # don't use IndexedComponent here, variables are always indexed components 
        # could also do this by iterating over component_objects(SimpleVar)...?
        if not isinstance(var,SimpleVar):
            for idx in var:
                if not (var[idx].lb is None):
                    if (var[idx].value < var[idx].lb - 1.0e-8):
                        print(var.name,idx)
                if not (var[idx].ub is None):
                    if (var[idx].value > var[idx].ub + 1.0e-8):
                        print(var.name,idx)
        else:
            if var.has_lb():
                if (var.value > var.ub + 1.0e-8):
                    print(var.name)
            if var.has_ub():
                if (var.value < var.lb - 1.0e-8):
                    print(var.name)
    print('- - -\n')

def write_differential_equations(flowsheet,suffix=''):

    m = flowsheet.MB_fuel

    with open('dCgdt_eqn'+suffix+'.txt','w') as f:
        m.eq_b1.pprint(ostream=f)
    with open('dqdt_eqn'+suffix+'.txt','w') as f:
        m.eq_b2.pprint(ostream=f)
    with open('dTgdt_eqn'+suffix+'.txt','w') as f:
        m.eq_d1.pprint(ostream=f)
    with open('dTsdt_eqn'+suffix+'.txt','w') as f:
        m.eq_d6.pprint(ostream=f)

    print('Time-differential equations written to files')

def get_vector_from_flowsheet(nlp,flowsheet):
    
    x = nlp.create_vector_x()
    order = nlp.variable_order()

    for i in range(0,nlp.nx):

        # extract "name" prefix from order (everything before '[')
        name = ''
        index_string = ''
        j = 0
        for char in order[i]:
            if char == '[' or char ==  '':
                j = j + 1
                # then create string containing the indices
                while order[i][j] != ']' and order[i][j] != '':
                    index_string = index_string + order[i][j] 
                    j = j + 1
                break
            else:
                name = name + char
            j = j + 1

        # create list of indices in the correct data types
        indices = []
        temp = ''
        for char in index_string:
            if char != ',' and char != '':
                temp = temp + char
            else:
                indices.append(temp)
                temp = ''
        if temp != '':
            indices.append(temp)

        for j in range(0,len(indices)):
            if indices[j] == '':
                raise ValueError('Did not expect this.')
            if indices[j][0].isdigit():
                indices[j] = float(indices[j])
        
        # evaluate the named variables at the correct indices
        if indices == []:
            x[i] = eval('flowsheet.'+name).value
        else:
            x[i] = eval('flowsheet.'+name+str(indices)).value

    return x
       
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

    write_differential_equations(flowsheet)

    # Then perturb
    solid_x_ptb = {'Fe2O3':0.25, 'Fe3O4':0.01, 'Al2O3':0.74}
    gas_y_ptb = {'CO2':0.03999, 'H2O':0.00001, 'CH4':0.96}
    #perturbInputs(flowsheet,0,Solid_M=691.4,Solid_T=1283,Solid_x=solid_x_ptb,
    #        Gas_F=150,Gas_T=350,Gas_y=gas_y_ptb)
    for t in mb.t:
        perturbInputs(flowsheet,t,Solid_M=691.4)

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
    mb.Tw_GW.fix(0.0)
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

    # other tentatively unused variables:
    mb.mFe_mAl.fix(0.0)
    mb.Solid_Out_M_Comp.fix()

    mb.eq_c5.deactivate()
    

    # Create a solver
    tol = 1e-8
    opt = SolverFactory('ipopt')
    opt.options = {'tol': tol,
                   'linear_solver' : 'ma57',
                   'bound_push': 1e-8,
                   'max_cpu_time': 600,
                   'print_level': 5}
                   #'halt_on_ampl_error': 'yes'}

    # initialized at steady state, works regardless:
    flowsheet.strip_bounds()

    #for z in mb.z:
    #    for t in mb.t:
    #        mb.Cg[z,'CH4',t].setlb(1e-8)

    for t in mb.t:
        alg_update(flowsheet,t)
        update_time_derivatives(flowsheet,t)


    print_violated_constraints(flowsheet,tol)

    nlp_ss = PyomoNLP(ss_flowsheet)
    x_ss = get_vector_from_flowsheet(nlp_ss,ss_flowsheet)
    jac_ss = nlp_ss.jacobian_g(x_ss)

    print('calculating steady state condition number...')
    ss_condition = np.linalg.cond(jac_ss.toarray())
    print('steady state condition number: ',ss_condition)

    fig1,ax1 = plt.subplots()
    ax1.jac_ss = plt.spy(jac_ss)
    ax1.set_facecolor('none')
    fig1.savefig('jac_ss.png',facecolor='none',edgecolor='none')#'#f2f2f2',edgecolor='none')



    nlp = PyomoNLP(flowsheet)
    v_order = nlp.variable_order()
    c_order = nlp.constraint_order()
    x = get_vector_from_flowsheet(nlp,flowsheet)
    lam = nlp.create_vector_y()

    jac_c = nlp.jacobian_g(x)
    hess_lag = nlp.hessian_lag(x,lam)
    kkt = BlockSymMatrix(2)
    kkt[0,0] = hess_lag
    kkt[1,0] = jac_c

    fig2,ax2 = plt.subplots()
    ax2.jac_c = plt.spy(jac_c)
    ax2.set_facecolor('none')
    fig2.savefig('jac_c.png',facecolor='none',edgecolor='none')


    #MA27 = hsl.MA27_LinearSolver()
    #jac_row_fortran = np.zeros(jac_c.nnz,dtype=np.intc) 
    #jac_col_fortran = np.zeros(jac_c.nnz,dtype=np.intc)
    #values = jac_c.data
    #for i in range(0,jac_c.nnz):
    #    jac_row_fortran[i] = int(jac_c.row[i] + 1)
    #    jac_col_fortran[i] = int(jac_c.col[i] + 1)
    #print('Doing symbolic factorization...')
    #MA27.DoSymbolicFactorization(nlp.nx,jac_row_fortran,jac_col_fortran)

    #print(jac_row_fortran)
    #print(jac_col_fortran)
    #print('Doing numeric factorization...')
    #num_status = MA27.DoNumericFactorization(nlp.nx,values)
    #print('Status: ',num_status)

    #jac_indices = range(0,jac_c.nnz)
    #for i in range(0,jac_c.nnz):
    #    if np.abs(values[i]) <= 1e-6:
    #        print('%0.2e'%values[i],str(jac_indices[i])+'-th nonzero.',jac_c.row[i],jac_c.col[i],
    #                c_order[jac_c.row[i]],v_order[jac_c.col[i]])

    #plot_switch = 0
    #if plot_switch == 1:
    #    fig,ax = plt.subplots()
    #    jac_value_plot = ax.bar(jac_indices,values)
    #    ax.set_yscale('log')
    #    fig.savefig('plots/jac_values.png')

    #print('calculating condition number...')
    #condition = np.linalg.cond(jac_c.toarray())
    #print('condition number: ',condition)

    #mb.input_objective = Objective(expr=sum((mb.Solid_In_M[t] -601.4)**2 for t in mb.t))

    flowsheet.write('fs_dyn.nl')

    #with open('factorized_fs.txt','w') as f:
    #    flowsheet.display(ostream=f)

    return flowsheet
    
if __name__ == "__main__":
    main()
