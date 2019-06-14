from __future__ import division
from __future__ import print_function

from pyomo.environ import value, Var, Constraint, Objective
from pyomo.environ import ConcreteModel
from pyomo.environ import Block
from pyomo.core.base.sets import _SetProduct
from pyomo.core.base.constraint import SimpleConstraint
from pyomo.core.base.var import SimpleVar
from pyomo.core.base.indexed_component import IndexedComponent
from pyomo.opt import SolverFactory
from pyomo.dae import DerivativeVar, ContinuousSet, Simulator

import time
import scipy
import casadi

from idaes_models.core import FlowsheetModel, ProcBlock
import MB_CLC as MB_CLC_fuel

import pdb

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
        fs.MB_fuel.Solid_In_Ts[t].fix(1183.15)      # kK
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
                    # this is the case of a variable indexed only by time
                    for t in time:
                        var_ol[t].set_value(ss_value)
                else:
                    var_ol.set_value(ss_value)
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

    print('Variable bounds violated')
    for var in flowsheet.MB_fuel.component_objects(Var,active=True):
        # don't use IndexedComponent here, variables are always indexed components 
        # could also do this by iterating over component_objects(SimpleVar)...?
        if not isinstance(var,SimpleVar):
            for idx in var:
                if not (var[idx].lb is None):
                    if (var[idx].value < var[idx].lb - 1.0e-8):
                        pdb.set_trace()
                        print(var.name,idx)
                if not (var[idx].ub is None):
                    if (var[idx].value > var[idx].ub + 1.0e-8):
                        pdb.set_trace()
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

def make_square(m):
    # should put this in a dedicated ~intialize~ function
    # that also intelligently initializes the model after perturbation
    m.eq_d4.deactivate()
    m.eq_d5.deactivate()
    m.eq_d8.deactivate()
    m.eq_d9.deactivate()
    m.eq_d10.deactivate()
    m.eq_g7.deactivate()
    m.eq_g8.deactivate()
    m.eq_g10.deactivate()
    m.eq_g11.deactivate()
    m.eq_g12.deactivate()
    m.eq_g13.deactivate()
    m.eq_g14.deactivate()
    m.eq_g4.deactivate()
    m.eq_g5.deactivate()
    m.eq_g2.deactivate()
    m.Tg_GW.fix(0.0)
    m.Tw_GW.fix(0.0)
    m.Tg_refractory.fix(0.0)
    m.Tw_Wamb.fix()
    m.Tw.fix()
    m.Nuw.fix()
    m.Nu_ext.fix()
    m.hw.fix()
    m.hext.fix()
    m.hext2.fix()
    m.U.fix()
    m.Uw.fix()
    m.Pr_ext.fix()
    m.Ra.fix()
    m.Re.fix()
    ###

    # other tentatively unused variables:
    m.mFe_mAl.fix(0.0)
    m.Solid_Out_M_Comp.fix()

    # choose how to calculate certain algebraic variables:
    m.eq_c5.deactivate()

#    m.strip_bounds()

def make_flowsheet(**kwargs):
    if 'press_drop' in kwargs:
        p_drop = kwargs['press_drop']
    if 't_dae_method' in kwargs:
        method = kwargs['t_dae_method']
    if 'horizon' in kwargs:
        H = kwargs['horizon']
    if 'ncp_t' in kwargs:
        ncp = kwargs['ncp_t']

    @ProcBlock("Flowsheet")
    class _Flowsheet(FlowsheetModel):
        # this Flowsheet class should be used to create single-fe flowsheet 
        # models for use in integration
        # These will be loaded into a full-horizon Flowsheet model created in the dyn_sim script
        def __init__(self, *args, **kwargs):
            FlowsheetModel.__init__(self,*args,**kwargs)
    
        def build(self):
            nfe = 6
            fe_a = 1/4.0
            fe_b = 0.2
            fe_set = [0,0.004]
            for i in range(1,nfe+1):
                if i < nfe*fe_a:
                    fe_set.append(i*fe_b/(nfe*fe_a))
                elif i == nfe:
                    fe_set.append(1)
                else:
                    fe_set.append(fe_b + (i-nfe*fe_a)*(1-fe_b)/(nfe*(1-fe_a)))
    
            self.MB_fuel = MB_CLC_fuel.MB(
                    parent=self,
                    z_dae_method = 'OCLR',
                    #t_dae_method = 'OCLR',
                    t_dae_method = method,
                    press_drop = p_drop,
                    fe_set = fe_set,
                    ncp = 3,
                    horizon = H,
                    nfe_t = 1,
                    ncp_t = ncp)

    fs = Flowsheet(name='MB_Model')

    return fs

def load_fe(fe,fs,t0):
    # t is the initial time of the single-finite-element 
    ms = fe.MB_fuel
    mt = fs.MB_fuel
    time = fe.MB_fuel.t
    time_fs = fs.MB_fuel.t
    for vfe in fe.MB_fuel.component_objects(Var,active=True):
        var_name = vfe.getname()
        vfs = getattr(fs.MB_fuel,var_name)

        if isinstance(vfe,SimpleVar):
            continue

        elif vfe.index_set() == time:
            for t in time:
                if t != time.first() or t0 == time_fs.first():
                    vfs[t+t0].set_value(vfe[t].value)

        elif vfe.index_set().dimen >= 2:
            if time not in vfe.index_set().set_tuple:
                continue
            else:
                for index in vfe:
                    t = index[vfe.index_set().dimen-1]

                    if t != time.first() or t0 == time_fs.first():
                        index_list = []

                        for i in range(0,vfe.index_set().dimen-1):
                            index_list.append(index[i])
                            # ^ append all non-time indices to the list
                        index_list.append(t+t0)

                        fs_index = tuple(index_list)
                        vfs[fs_index].set_value(vfe[index].value)


