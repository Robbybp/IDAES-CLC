from __future__ import division
from __future__ import print_function
from pyomo.environ import *
from pyomo.core.base.constraint import SimpleConstraint
from pyomo.core.base.var import SimpleVar, IndexedVar
from pyomo.dae import *
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes_models.core import FlowsheetModel, ProcBlock
import pdb
import ss_sim 
from idaes_models.core import FlowsheetModel, ProcBlock
import mb_clc as MB_CLC_fuel
from utils import setInputs, perturbInputs, setICs, initialize_ss, make_flowsheet, make_square, print_violated_constraints
#import dyn_sim 

# ^ for debugging purposes 

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
                t_dae_method = 'BFD1',
                press_drop = 'Ergun',
                fe_set = fe_set,
                ncp = 3,
                horizon = 5,
                nfe_t = 1,
                ncp_t = 1)


def void_alg_variables(fs,t):
    m = fs.MB_fuel

    ###
    ## shouldn't do this because it will be applied for every point in time...
    ## ... could fix to void only the algebraic variables 
    ## should do this before performing an algebraic update for all time... 
    dif_var_list = [ m.Cg.name, m.q.name, m.Tg.name, m.Ts.name ]
    for var in m.component_objects(Var):
        
        # want some logic here like 'if time-index is t:'
        # or even better, if time-index is not t: continue'
        # first step: is variable indexed by time?
        fixed = True 
        for index in var:
            # if any of the coordinates of the variable are unfixed, then consider the variable unfixed
            # otherwise it is fixed, and disregard it as a possible algebraic variable
            if var[index].fixed == False:
                fixed = False
                break

        stale = True
        for index in var: 
            if var[index].stale == False:
                stale = False
                break
        
        is_time_derivative = False
        # this checks if its a derivative variable indexed by time ?
        # or if its a derivative variable with respect to time?
        if isinstance(var,DerivativeVar):
            if m.t in var.get_continuousset_list():
                is_time_derivative = True

        if fixed == False and (var.name not in dif_var_list) and is_time_derivative == False:
            if isinstance(var,SimpleVar):
                continue
            # if not a simple var, then index_set is defined
            elif var.index_set() == m.t:
                var[t].set_value(None)
            elif var.index_set().dimen > 1:
                if m.t in var.index_set().set_tuple:
                    for index in var:
                        n = var.index_set().dimen
                        if index[n-1] == t:
                            var[index].set_value(None)

def alg_update(fs,t):
    # function to update algebraic variables of the model at a specified time step, placeholder for now
    # this should be the first thing to do, because everything else I do here will rely on it
    m = fs.MB_fuel

    ###
    ## shouldn't do this because it will be applied for every point in time...
    ## ... could fix to void only the algebraic variables 
    ## should do this before performing an algebraic update for all time... 
    #dif_var_list = [ m.Cg.name, m.q.name, m.Tg.name, m.Ts.name ]
    #for var in m.component_objects(Var):
    #    
    #    # want some logic here like 'if time-index is t:'
    #    # or even better, if time-index is not t: continue'
    #    # first step: is variable indexed by time?
    #    fixed = True 
    #    for index in var:
    #        # if any of the coordinates of the variable are unfixed, then consider the variable unfixed
    #        # otherwise it is fixed, and disregard it as a possible algebraic variable
    #        if var[index].fixed == False:
    #            fixed = False
    #            break

    #    stale = True
    #    for index in var: 
    #        if var[index].stale == False:
    #            stale = False
    #            break
    #    
    #    is_time_derivative = False
    #    # this checks if its a derivative variable indexed by time ?
    #    # or if its a derivative variable with respect to time?
    #    if isinstance(var,DerivativeVar):
    #        if m.t in var.get_continuousset_list():
    #            is_time_derivative = True

    #    if fixed == False and (var.name not in dif_var_list) and is_time_derivative == False:
    #        if isinstance(var,SimpleVar):
    #            continue
    #        # if not a simple var, then index_set is defined
    #        elif var.index_set() == m.t:
    #            var[t].set_value(None)
    #        elif var.index_set().dimen > 1:
    #            if m.t in var.index_set().set_tuple:
    #                for index in var:
    #                    n = var.index_set().dimen
    #                    if index[n-1] == t:
    #                        var[index].set_value(None)

            #else:
            #    for index in var:
            #        var[index].set_value(None)
    ###
    
    #
    # "differential" variables that are specified by (inlet) boundary conditions
    # (along with solid velocity)
    #

    calculate_variable_from_constraint(m.Tg[0,t],m.eq_f3[t])
    calculate_variable_from_constraint(m.Ts[1,t],m.eq_f4[t])

    # vs:
    calculate_variable_from_constraint(m.vs[t],m.eq_f6[t])
    # vg[0]:
    m.MW_vap[0,t].set_value(1e-3*sum( m.Gas_In_y[j,t].value*m.MW[j] for j in m.GasList ) )
    m.P[0,t].set_value(m.Gas_In_P[t].value)
    calculate_variable_from_constraint(m.rho_vap[0,t],m.eq_p6[0,t])
    calculate_variable_from_constraint(m.vg_in[t],m.eq_f5[t])
    m.vg[0,t].set_value(m.vg_in[t].value)

    for j in m.GasList:
        calculate_variable_from_constraint(m.G_flux[0,j,t],m.eq_f1[j,t])
        calculate_variable_from_constraint(m.F[0,j,t],m.eq_c1[0,j,t])
        calculate_variable_from_constraint(m.Cg[0,j,t],m.eq_c4[0,j,t])
    for j in m.SolidList:
        calculate_variable_from_constraint(m.S_flux[1,j,t],m.eq_f2[j,t])
        calculate_variable_from_constraint(m.q[1,j,t],m.eq_c13[1,j,t])

    #
    # calculating algebraic variables from constraints, with differential variables and inputs/disturbances known 
    #

    for z in m.z:
        # CgT:
        m.CgT[z,t].set_value( sum( m.Cg[z,j,t].value for j in m.GasList ) )
        # y:
        for j in m.GasList:
            m.y[z,j,t].set_value( m.Cg[z,j,t].value/m.CgT[z,t].value )
        # ytot:
        m.ytot[z,t].set_value( sum(m.y[z,j,t].value for j in m.GasList) )
        # mu_vap:
        calculate_variable_from_constraint(m.mu_vap[z,t],m.eq_p8[z,t])
        # MW_vap:
        calculate_variable_from_constraint(m.MW_vap[z,t],m.eq_p5[z,t])

        #if z == 0:
        #    m.P[z,t].set_value(m.Gas_In_P[t].value)
        #else:
        #    m.P[z,t].set_value( m.CgT[z,t].value * m.R.value*1e-5 * m.Tg[z,t].value )

        # P:
        calculate_variable_from_constraint(m.P[z,t],m.eq_e1[z,t])
        # rho_vap:
        calculate_variable_from_constraint(m.rho_vap[z,t],m.eq_p6[z,t])

    # vg_in:
    calculate_variable_from_constraint(m.vg_in[t],m.eq_f5[t])
    # could do this one finite element at a time, directly after initializing P-values
    verbose_vg = False 
    for z in m.z:    
        if z != 0:
            prev = 0
            for temp in m.z:
                if temp < z:
                    prev = temp
                else:
                    break
            v_prev = (m.vg[prev,t].value+(1e-3*m.vs[t].value))

            # dPdz:
            calculate_variable_from_constraint(m.dPdz[z,t],m.dPdz_disc_eq[z,t])
            # vg:
            if verbose_vg: print('\t',z)
            if verbose_vg: print('\tvg pre-update:\t',m.vg[z,t].value)
            if verbose_vg: print('\tdPdz: ',m.dPdz[z,t].value)
            a = 1.75/150*m.rho_vap[z,t].value*(1-m.eps.value)*m.L.value/\
                    (m.dp.value*m.eps.value**3)
            b = (1e-3*m.mu_vap[z,t].value)*(1-m.eps.value)**2*m.L.value/\
                    (m.dp.value**2*m.eps.value**3)
            c = 1e5*m.dPdz[z,t].value/150

            d = (b**2-4*a*c)

            if a < 0:
                raise ValueError('Sanity check failed, a < 0 initially')
            vg_negative = False
            if d < 0:
                a = -a
                d = (b**2-4*a*c)
            if d < 0:
                # if discriminant is still negative after changing sign,
                raise ValueError('Cannot solve for vg - quadratic has no real root')
            vp = (-b + sqrt(d))/(2*a)
            vm = (-b - sqrt(d))/(2*a)
            if a*vp < 0 and a*vm < 0:
                raise ValueError('(vs+vg) does not have consistent sign')
            elif a*vp < 0:
                v = vm
            elif a*vm < 0:
                v = vp
            elif abs(vm-v_prev) < abs(vp-v_prev):
                v = vm
            else:
                v = vp

            m.vg[z,t].set_value(v-(1e-3*m.vs[t].value))
            # if calculation fails, should try switching sign of v_g+v_s
            #calculate_variable_from_constraint(m.vg[z,t],m.eq_e2[z,t])#,linesearch=False)
            if verbose_vg: print('\tvg post-update:\t',m.vg[z,t].value)
        
    # Gas_Out_P:
    calculate_variable_from_constraint(m.Gas_Out_P[t],m.eq_f7[t])

    # S_flux:
    # Solid_M:
    # Solid_M_total:
    # Solid_Out_M:
    #calculate_variable_from_constraint(m.Solid_Out_M[t],m.eq_f8[t])
    # Solid_Out_M_Comp is stale... does it need a value?
    #calculate_variable_from_constraint(m.Solid_Out_M_Comp

    # Solid_Out_Ts:
    #calculate_variable_from_constraint(m.Solid_Out_Ts[t],m.eq_f9[t])
    # Solid_Out_x:
    #for j in m.SolidList:
    #    calculate_variable_from_constraint(m.Solid_Out_x[j,t],m.eq_f10[j,t])

    # Pr_ext:
    #calculate_variable_from_constraint(m.Pr_ext,m.eq_g4)
    # ^ if it's not indexed by time, I don't care about it

    for z in m.z:
        # qT:
        m.qT[z,t].set_value( sum( m.q[z,j,t].value for j in m.SolidList ) )
        # x:
        for j in m.SolidList:
            m.x[z,j,t].set_value( m.q[z,j,t].value/m.qT[z,t].value )
        # xtot:
        m.xtot[z,t].set_value( sum(m.x[z,j,t].value for j in m.SolidList) )
        # X:
        calculate_variable_from_constraint(m.X[z,t],m.eq_r2[z,t])
        if z != 0:
            # k:
            for i in m.rxn_idx:
                calculate_variable_from_constraint(m.k[z,i,t],m.eq_r1[z,i,t])
            # X_term:
            m.X_term[z,t].set_value( (1-m.X[z,t].value)**(2/3) )
            calculate_variable_from_constraint(m.X_term[z,t],m.eq_r3[z,t])
            # r_gen:
            for i in m.rxn_idx:
                calculate_variable_from_constraint(m.r_gen[z,i,t],m.eq_r4[z,i,t])
            # rs, rg:
            for j in m.SolidList:
                calculate_variable_from_constraint(m.rs[z,j,t],m.eq_r5[z,j,t])
            for j in m.GasList:
                calculate_variable_from_constraint(m.rg[z,j,t],m.eq_r6[z,j,t])
            # Ctrans, qtrans:
            for j in m.GasList:
                calculate_variable_from_constraint(m.Ctrans[z,j,t],m.eq_b3[z,j,t])
            for j in m.SolidList:
                calculate_variable_from_constraint(m.qtrans[z,j,t],m.eq_b4[z,j,t])

        # Gas flows:
        # F:
        if z != m.z.first():
            for j in m.GasList:
                calculate_variable_from_constraint(m.F[z,j,t],m.eq_c4[z,j,t])
        # Ftotal:
        m.Ftotal[z,t].set_value( sum( m.F[z,j,t].value for j in m.GasList) )
        # G_flux:
        if z != m.z.first():
            for j in m.GasList:
                calculate_variable_from_constraint(m.G_flux[z,j,t],m.eq_c1[z,j,t])
        #else:
        #    for j in m.GasList:
        #        calculate_variable_from_constraint(m.F[z,j,t],m.eq_c1[z,j,t])

        # Gas_M:
        for j in m.GasList:
            calculate_variable_from_constraint(m.Gas_M[z,j,t],m.eq_c2[z,j,t])

        # Solid flows:
        # S_flux:
        if z != m.z.last():
            for j in m.SolidList:
                calculate_variable_from_constraint(m.S_flux[z,j,t],m.eq_c13[z,j,t])
        # qT and x, already done
        # Solid_M:
        for j in m.SolidList:
            calculate_variable_from_constraint(m.Solid_M[z,j,t],m.eq_c9[z,j,t])
        # Solid_M_total:
        m.Solid_M_total[z,t].set_value( 
                sum(m.Solid_M[z,j,t].value for j in m.SolidList) )
        # Solid_F:
        for j in m.SolidList:
            calculate_variable_from_constraint(m.Solid_F[z,j,t],m.eq_c11[z,j,t])
        m.Solid_F_total[z,t].set_value(
                sum(m.Solid_F[z,j,t].value for j in m.SolidList) )

        # Temperatures and stuff:

        # k_vap:
        calculate_variable_from_constraint(m.k_vap[z,t],m.eq_p16[z,t])
        # cp_vap:
        calculate_variable_from_constraint(m.cp_vap[z,t],m.eq_p10[z,t])
        # cp_gas:
        calculate_variable_from_constraint(m.cp_gas[z,t],m.eq_p11[z,t])
        # Pr:
        calculate_variable_from_constraint(m.Pr[z,t],m.eq_g3[z,t])
        # Rep:
        calculate_variable_from_constraint(m.Rep[z,t],m.eq_g1[z,t])
        # Nu:
        calculate_variable_from_constraint(m.Nu[z,t],m.eq_g6[z,t])
        # hf:
        calculate_variable_from_constraint(m.hf[z,t],m.eq_g9[z,t])
        # Tg_GS:
        calculate_variable_from_constraint(m.Tg_GS[z,t],m.eq_d3[z,t])

        # DH_rxn_s, Ts_dHr:
        if z != m.z.first():
            calculate_variable_from_constraint(m.DH_rxn_s[z,t],m.eq_p2[z,t])
            calculate_variable_from_constraint(m.Ts_dHr[z,t],m.eq_d7[z,t])

        ### variables only used for calculation of eternal fluxes need not be calculated, 
        ### as these fluxes are fixed to zero for now. 
        ### if external fluxes are to be used, need to solve a system to 
        ### determine the mutually dependent Tw[t,1] and Ra[t]
        # Nuw:
        #calculate_variable_from_constraint(m.Nuw[z,t],m.eq_g7[z,t])
        # hw:
        #calculate_variable_from_constraint(m.hw[z,t],m.eq_g10[z,t])
        # Ra:
        # Nu_ext:
        # 
        # Tw_GW:
        ###

        # umf: 
        calculate_variable_from_constraint(m.umf[z,t],m.eq_e3[z,t])
        # v_diff:
        calculate_variable_from_constraint(m.v_diff[z,t],m.eq_e4[z,t])
        # Gh_flux:
        calculate_variable_from_constraint(m.Gh_flux[z,t],m.eq_d2[z,t])
        # cp_sol:
        calculate_variable_from_constraint(m.cp_sol[z,t],m.eq_p4[z,t])
        # Sh_flux:
        calculate_variable_from_constraint(m.Sh_flux[z,t],m.eq_d11[z,t])
        # cv_vap:
        calculate_variable_from_constraint(m.cv_vap[z,t],m.eq_p12[z,t])
        # k_cpcv:
        calculate_variable_from_constraint(m.k_cpcv[z,t],m.eq_p13[z,t])

    # derivatives
    for z in m.z:
        if z != m.z.first():
            for j in m.GasList:
                calculate_variable_from_constraint(m.dG_fluxdz[z,j,t],m.dG_fluxdz_disc_eq[z,j,t])
            for j in m.SolidList:
                calculate_variable_from_constraint(m.dS_fluxdz[z,j,t],m.dS_fluxdz_disc_eq[z,j,t])
            calculate_variable_from_constraint(m.dGh_fluxdz[z,t],m.dGh_fluxdz_disc_eq[z,t])
            calculate_variable_from_constraint(m.dSh_fluxdz[z,t],m.dSh_fluxdz_disc_eq[z,t])
            #calculate_variable_from_constraint(m.dTsdz[z,t],m.dTsdz_disc_eq[z,t])
        
    # outlet and overall parameters, not dependent on z:

    # X_gas:
    calculate_variable_from_constraint(m.X_gas[t],m.eq_c8[t])
    # X_OC:
    calculate_variable_from_constraint(m.X_OC[t],m.eq_c17[t])
    # Solid_Out_M:
    calculate_variable_from_constraint(m.Solid_Out_M[t],m.eq_f8[t])
    # Solid_Out_M_Comp: (stale, literally does not get used, 
    # no equation to specify from ...)
    #calculate_variable_from_constraint(m.Solid_Out_M_Comp...
    # Solid_Out_Ts:
    calculate_variable_from_constraint(m.Solid_Out_Ts[t],m.eq_f9[t])
    # Solid_Out_x:
    for j in m.SolidList:
        calculate_variable_from_constraint(m.Solid_Out_x[j,t],m.eq_f10[j,t])
    
    #with open('post_alg_update.txt','w') as f:
    #    m.display(ostream=f)


def update_time_derivatives(fs,t):
    verbose = False 
    debug = True
    m = fs.MB_fuel
    if verbose: print('\nTime derivative values: \n - - -')

    for z in m.z:
        if z != m.z.first():
            for j in m.GasList:
                #if t > 0: pdb.set_trace()
                calculate_variable_from_constraint(m.dCgdt[z,j,t],m.eq_b1[z,j,t])
                if verbose: print('\t\tdCgdt ',z,',',j,',',t,':\t',m.dCgdt[z,j,t].value)

            for j in m.SolidList:
                calculate_variable_from_constraint(m.dqdt[z,j,t],m.eq_b2[z,j,t])
                if verbose: print('\t\tdqdt ',z,',',j,',',t,':\t',m.dqdt[z,j,t].value)

            calculate_variable_from_constraint(m.dTgdt[z,t],m.eq_d1[z,t])
            if verbose: print('\t\tdTgdt ',z,',',t,':\t',m.dTgdt[z,t].value)

            calculate_variable_from_constraint(m.dTsdt[z,t],m.eq_d6[z,t])
            if verbose: print('\t\tdTsdt ',z,',',t,':\t',m.dTsdt[z,t].value)

        #else:
        #    for j in m.GasList:
        #        m.dCgdt[z,j,t].set_value(0)
        #    for j in m.SolidList:
        #        m.dqdt[z,j,t].set_value(0)
        #    m.dTgdt[z,t].set_value(0)
        #    m.dTsdt[z,t].set_value(0)
    if verbose: print('- - -\n')

def ee_update(fs,t):
    # performs explicit euler update (starting at t)
    # to the differential variables 

    # assumes the values of the differential variables and their time derivatives are correct at t

    verbose = False 
    residual = True
    m = fs.MB_fuel
    tp = t
    for ti in m.t:
        if ti > t:
            tp = ti
            break
    if t == tp:
        raise ValueError('Cannot perform exp. Euler update to last time point')

    h = tp - t
    dif_var_list = [ m.Cg, m.q, m.Tg, m.Ts]
    if verbose: print('differential variables: \n- - -')
    for var in dif_var_list:
        dvar_name = 'm.' + 'd' + var.local_name + 'dt'
        dvar = eval(dvar_name)

        if var.index_set() == m.t:
            # if the variable is only indexed by time the euler update
            # only needs to be performed once
            var[tp].set_value( var[t].value + h*dvar[t].value )
            if verbose: print('\t\t',var.name,':\t',var[t].value,' ---> ',var[tp].value)

        else:
            # the euler update needs to be performed for each non-time
            # index of the variable
            n = var.index_set().dimen
            sets = var.index_set().set_tuple
            # temp set will contain the non-time indices of the variable
            # note: set_tuple has at least two elements (__ and m.t)
            temp_set = sets[0]

            if n == 2:
                # temp_set is the non-time indexing set for var
                for index in temp_set:
                    time = (index,t)
                    plus = (index,tp)

                    new_value = var[time].value + h*dvar[time].value
                    if var[plus].has_lb():
                        if new_value < var[plus].lb:
                            new_value = var[plus].lb 
                    var[plus].set_value( new_value )

                    if verbose: print('\t\t',var.name,':\t',var[time].value,' ---> ',var[plus].value)
                    if verbose and residual: print('\t\t\t\t',var[plus].value-var[time].value)

            if n > 2:
                # here temp_set is the product of all non-time indexing sets
                for i in range(1,n-1):
                    temp_set = temp_set * sets[i]
                for index in temp_set:
                    # indices are already tuples
                    time = index + (t,)
                    plus = index + (tp,)

                    new_value = var[time].value + h*dvar[time].value
                    if var[plus].has_lb():
                        if new_value < var[plus].lb:
                            new_value = var[plus].lb
                    var[plus].set_value( new_value )

                    if verbose: print('\t\t',var.name,':\t',var[time].value,' ---> ',var[plus].value)
                    if verbose and residual: print('\t\t\t\t',var[plus].value-var[time].value)

def propagate_alg_variables(fs,t):
    # initialize the values of alg vars at t+ with those at t,
    # this is an initialization for the (maybe implicity) update_alg function
    m = fs.MB_fuel
    tp = t
    for ti in m.t:
        if ti > t:
            tp = ti
            break
    if t == tp:
        raise ValueError('Cannot propagate algebraic variables beyond last time point')
    # find the algebraic variables (that are indexed by time)
    # separate into indexed by only time, indexed by time and one other set,
    #     and indexed by time and 2+ other sets
    # var[(...,t+)].set_value(var[(...,t)])
                
    dif_var_list = [ m.Cg.name, m.q.name, m.Tg.name, m.Ts.name ]
    for var in m.component_objects(Var):
        
        fixed = True 
        for index in var:
            # if any of the coordinates of the variable are unfixed, then consider the variable unfixed
            # otherwise it is fixed, and disregard it as a possible algebraic variable
            if var[index].fixed == False:
                fixed = False
                break

        stale = True
        for index in var: 
            if var[index].stale == False:
                stale = False
                break
        
        is_time_derivative = False
        # this checks if its a derivative variable indexed by time ?
        # or if its a derivative variable with respect to time?
        if isinstance(var,DerivativeVar):
            if m.t in var.get_continuousset_list():
                is_time_derivative = True

        if fixed == False and (var.name not in dif_var_list) and is_time_derivative == False:
            if isinstance(var,SimpleVar):
                continue
            # if not a simple var, then index_set is defined
            elif var.index_set() == m.t:
                var[tp].set_value( var[t].value )
            elif var.index_set().dimen == 1:
                # if only indexed by one set, and it's not time, I don't care about it 
                continue
            # now know that var is indexed by more than one set, so set_tuple is defined
            elif m.t not in var.index_set().set_tuple:
                # if not indexed by time, I don't care about it
                continue
            elif var.index_set().dimen == 2:
                # indexed by two sets, the second of which (by convention, here) is time
                sets = var.index_set().set_tuple
                temp_index_set = sets[0]
                #print(var,sets)
                for i in temp_index_set:
                    var[(i,tp)].set_value( var[(i,t)].value )
            elif var.index_set().dimen > 2:
                sets = var.index_set().set_tuple
                temp_index_set = sets[0]
                n = var.index_set().dimen
                #print(var,sets)
                for i in range(1,n-1):
                    temp_index_set = temp_index_set * sets[i]
                for index in temp_index_set:
                    # indices are already tuples
                    time = index + (t,)
                    plus = index + (tp,)
                    var[plus].set_value( var[time].value )

def initialize_next_t(fs):
    m = fs.MB_fuel

    # logic:
    # for all variables in model:
    #   if indexed by time:
    #       set_value( variable[t.last()] )
    # (really only need to propagate differential variables/inputs,
    # as algebraic variables and time derivatives will be updated)
    for var in m.component_objects(Var):
        if isinstance(var,SimpleVar):
            continue

        elif var.index_set == m.t:
            for t in m.t:
                last = m.t.last()
                if t != last:
                    var[t].set_value(last)

        elif var.index_set().dimen > 1:
            sets = var.index_set().set_tuple

            if m.t not in sets:
                continue

            n = var.index_set().dimen
            if n == 2 and m.t == sets[n-1]:
                partial_index = sets[0]
                for index in partial_index:
                    next_fe = (index,m.t.last())
                    for t in m.t:
                        if t != m.t.last():
                            var[(index,t)].set_value( var[next_fe].value )
            elif n > 2 and m.t == sets[n-1]:
                partial_index = sets[0]
                for i in range(1,n-1):
                    partial_index = partial_index * sets[i]
                for index in partial_index:
                    next_fe = index + (m.t.last(),)
                    for t in m.t:
                        if t != m.t.last():
                            var[index + (t,)].set_value( var[next_fe].value )


def implicit_integrate(fs):
    # expects fs to be a model defined over a single finite element. 
    # add option to this function for initialization method (?)
    m = fs.MB_fuel

    for t in m.t:
        alg_update(fs,t)
        update_time_derivatives(fs,t)

    # initialize, at first via explicit integration 

    opt = SolverFactory('ipopt')    
    opt.options = {'tol': 1e-8,
                   'linear_solver' : 'ma57',
                   'bound_push': 1e-8,
                   'max_cpu_time': 600,
                   'print_level': 5,
                   'linear_system_scaling': 'none',
                   #'linear_scaling_on_demand': 'no',
                   #'nlp_scaling_method': 'none',
                   #'nlp_scaling_max_gradient': 1e6,
                   'halt_on_ampl_error': 'yes'}
    
    results = opt.solve(fs,tee=True)
    #fs_solved = fs.clone()
    #initialize_next_fe(fs)
    #return fs_solved

def ee_integrate(fs,t):
    alg_update(fs,t)
    update_time_derivatives(fs,t)
    ee_update(fs,t)
    propagate_alg_variables(fs,t)



    



def integrate(fs,t): 
    # function to integrate the model one time step, from t
    m = fs.MB_fuel

    # exp. euler integration:
    #
    # solve steady state, initialize t = 0
    # apply perturbation to t = 0
    # for collocation points:
    #   algebraic update
    #   update derivatives
    #   update differential variables at next time step
    #   (optional: apply new perturbation)

    # this function will:
    #        i)  apply algebraic update
    #       ii)  update derivatives
    #      iii)  update next set of differential variables
    #            (so should not be applied at last time point)

    # for an implicit method, need to solve for states at entire finite element simultaneously.
    # for explicit method, can solve for states at each time point individually
    # this function will do both. integrate explicitly, and use to initialize an
    # implicit integration.

    # void_alg_variables(fs,t) # or maybe have alg_update call void_alg_variables()
                               # but don't want to do this every time because it cancels out propagation...
    alg_update(fs,t)
    update_time_derivatives(fs,t)
    ee_update(fs,t)
    propagate_alg_variables(fs,t)

    # need to perform alg_update on fs at t.last() after calling integrate() 
    # at second-to-last time point, as it will not be called here


def set_fe_ICs(fe1,fe2):
    m1 = fe1.MB_fuel
    m2 = fe2.MB_fuel
    diff_vars = []
    diff_vars.append('Cg')
    diff_vars.append('q')
    diff_vars.append('Tg')
    diff_vars.append('Ts')

    for z in m1.z:
        for j in m1.GasList:
            m2.Cg_0[z,j].set_value( m1.Cg[z,j, m1.t.last() ].value )
        for j in m1.SolidList:
            m2.q_0[z,j].set_value( m1.q[z,j, m1.t.last() ].value )
        m2.Tg_0[z].set_value( m1.Tg[z, m1.t.last() ].value )
        m2.Ts_0[z].set_value( m1.Ts[z, m1.t.last() ].value )

def initialize_next_fe(fe1,fe2):
    
    time1 = fe1.MB_fuel.t
    time2 = fe2.MB_fuel.t
    for var1 in fe1.MB_fuel.component_objects(Var,active=True):
        var_name = var1.getname()
        var2 = getattr(fe2.MB_fuel,var_name)

        if isinstance(var1,SimpleVar):
            var2.set_value(var1.value)
            continue


        if var2.index_set().dimen >= 2:
            if time2 not in var2.index_set().set_tuple:
                for index in var2:
                    var2[index].set_value( var1[index].value )
        elif var2.index_set() != time2:
            for index in var2:
                var2[index].set_value( var1[index].value )

        if var2.index_set() == time2:
            for t in time2:
                var2[t].set_value( var1[time1.last()].value )

        if var2.index_set().dimen >= 2:
            if time2 in var2.index_set().set_tuple:
                for index in var2:
                    index1_list = []
                    for i in range(0,var2.index_set().dimen - 1):
                        index1_list.append( index[i] )
                    index1_list.append( time1.last() )
                    index1 = tuple(index1_list)
                    var2[index].set_value( var1[index1].value )

def copy_inputs(mb, fs_list):
    fs_list[0].MB_fuel.Solid_In_M[0].fix(mb.Solid_In_M[0].value)
    fs_list[0].MB_fuel.Solid_In_Ts[0].fix(mb.Solid_In_Ts[0].value)
    for i in mb.SolidList:
        fs_list[0].MB_fuel.Solid_In_x[i,0].fix(mb.Solid_In_x[i,0].value)
    fs_list[0].MB_fuel.Gas_In_F[0].fix(mb.Gas_In_F[0].value)
    fs_list[0].MB_fuel.Gas_In_P[0].fix(mb.Gas_In_P[0].value)
    fs_list[0].MB_fuel.Gas_In_Tg[0].fix(mb.Gas_In_Tg[0].value)
    for i in mb.GasList:
        fs_list[0].MB_fuel.Gas_In_y[i,0].fix(mb.Gas_In_y[i,0].value)

    for fe in fs_list:
        for t in fs_list[fe].MB_fuel.t:
            if t != fs_list[fe].MB_fuel.t.first():
                fs_list[fe].MB_fuel.Solid_In_M[t].fix(mb.Solid_In_M[fe+t].value)
                fs_list[fe].MB_fuel.Solid_In_Ts[t].fix(mb.Solid_In_Ts[fe+t].value)
                for i in mb.SolidList:
                    fs_list[fe].MB_fuel.Solid_In_x[i,t].fix(mb.Solid_In_x[i,fe+t].value)
                fs_list[fe].MB_fuel.Gas_In_F[t].fix(mb.Gas_In_F[fe+t].value)
                fs_list[fe].MB_fuel.Gas_In_P[t].fix(mb.Gas_In_P[fe+t].value)
                fs_list[fe].MB_fuel.Gas_In_Tg[t].fix(mb.Gas_In_Tg[fe+t].value)
                for i in mb.GasList:
                    fs_list[fe].MB_fuel.Gas_In_y[i,t].fix(mb.Gas_In_y[i,fe+t].value)


def clc_integrate(mb):

    ss_flowsheet = ss_sim.main()

    time_set = mb.t.get_finite_elements()
    ncp = 0
    for t in mb.t:
        if t < time_set[1]:
            ncp = ncp + 1
        else:
            break
    if ncp == 1:
        method = 'BFD1'
    elif ncp > 1:
        method = 'OCLR'
    else: 
        raise ValueError('Bad value of ncp when parsing full-horizon model')

    #time_set = [0,5,10,15,20,25,30]
            #35,40,45,50,55,60,65,70,
            #75,80,85,90,95,100,105,110,115,120,
            #125,130,135,140,145,150,155,160,165]

    # Create a custom grid, fe_set        
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

    fs_list = {}
    for i in range(0,len(time_set)):
        if i != len(time_set)-1:
        # assumption here that the last element of the time_set
        # is the end of the horizon. Could this ever not be the case?
            H = time_set[i+1] - time_set[i]
            fs_list[time_set[i]] = make_flowsheet(press_drop ='Ergun',
                                        t_dae_method = method,
                                        horizon = H,
                                        ncp_t = ncp)
#            fs_list[t] = MB_CLC_fuel.MB(
#                parent = None,
#                z_dae_method = 'OCLR',
#                #t_dae_method = 'OCLR',
#                t_dae_method = method,
#                press_drop = 'Ergun',
#                fe_set = fe_set,
#                ncp = 3,
#                horizon = H,
#                nfe_t = 1,
#                ncp_t = ncp)
    
    time_set.pop()
    
    # set initial conditions to steady state
    for fe in time_set:
        setICs(fs_list[fe],ss_flowsheet)

    # initialize variables to steady state values
    for fe in time_set:
        initialize_ss(fs_list[fe],ss_flowsheet)

    # set input values for each finite element model
    for fe in time_set:
        setInputs(fs_list[fe])

    gas_y_ptb = {'CO2': 0.04999, 'H2O': 0.00001, 'CH4': 0.95}

    # perturb inputs
    copy_inputs(mb, fs_list)
    #for fe in time_set:
    #    for t in fs_list[fe].MB_fuel.t:
    #        copy_inputs(fs_list[fe].MB_fuel, t, mb, fe+t)
            #perturbInputs(fs_list[fe], t, Gas_y=gas_y_ptb)

    for t in time_set:
        make_square(fs_list[t].MB_fuel)
        

    #with open('one_fe.txt','w') as f:
    #    fs_list[0].display(ostream=f)

    for i in range(0,len(time_set)):
        # should initialization occur here or in the imp_int function?
        implicit_integrate(fs_list[time_set[i]])
        if i != len(time_set)-1:
            set_fe_ICs(fs_list[time_set[i]],fs_list[time_set[i+1]])
            initialize_next_fe(fs_list[time_set[i]],fs_list[time_set[i+1]])
            
            # set initial conditions of next fe to final values of current fe
            # initialize next fe, somehow
            # a) initialize to initial conditions b) continuation c) explicit integrator

    fe_width = time_set[0]
    # ^ assumes a fixed finite element width
    print([ fs_list[t].MB_fuel.Cg[0.00062,'H2O',fe_width].value for t in time_set ] )
        
    
    return fs_list

    

def main():
    # main routine
    print('executing main routine')

if __name__ == '__main__':
    main()
