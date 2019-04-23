from __future__ import division
from __future__ import print_function
from pyomo.environ import *
from pyomo.core.base.constraint import SimpleConstraint
from pyomo.core.base.var import SimpleVar
from pyomo.dae import *
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes_models.core import FlowsheetModel, ProcBlock
import MB_CLC as MB_CLC_fuel
import pdb

def alg_update(fs,t):
    # function to update algebraic variables of the model at a specified time step, placeholder for now
    # this should be the first thing to do, because everything else I do here will rely on it
    m = fs.MB_fuel
    dif_var_list = [ m.Cg.name, m.q.name, m.Tg.name, m.Ts.name ]
    #dif_var_list = [ m.Cg, m.q, m.Tg, m.Ts ]
    print(dif_var_list)
    for var in m.component_objects(Var):
        
        fixed = True 
        for index in var:
            # if any of the coordinates of the variable are unfixed, then consider the variable unfixed
            # otherwise it is fixed, and disregard it as a possible algebraic variable
            if var[index].fixed == False:
                fixed = False
                break
        
        is_time_derivative = False
        if isinstance(var,DerivativeVar):
            if m.t in var.get_continuousset_list():
                is_time_derivative = True

        if fixed == False and (var.name not in dif_var_list) and is_time_derivative == False:
            if isinstance(var,SimpleVar):
                var.set_value(None)
            else:
                for index in var:
                    var[index].set_value(None)
            print(var.local_name)

    # calculating algebraic variables from constraints, with differential variables and inputs/disturbances known 
    for z in m.z:
        # CgT:
        m.CgT[z,t].set_value( sum( m.Cg[z,j,t].value for j in m.GasList ) )
        # y:
        for j in m.GasList:
            m.y[z,j,t].set_value( m.Cg[z,j,t].value/m.CgT[z,t].value )
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

    # A_bed:
    calculate_variable_from_constraint(m.A_bed,m.eq_a2)
    # vg_in:
    calculate_variable_from_constraint(m.vg_in[t],m.eq_f5[t])
    # vs:
    calculate_variable_from_constraint(m.vs[t],m.eq_f6[t])
    # could do this one finite element at a time, directly after initializing P-values
    for z in m.z:    
        if z != 0:
            # dPdz:
            calculate_variable_from_constraint(m.dPdz[z,t],m.dPdz_disc_eq[z,t])
        # vg:
        calculate_variable_from_constraint(m.vg[z,t],m.eq_e2[z,t])
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

    for z in m.z:
        m.qT[z,t].set_value( sum( m.q[z,j,t].value for j in m.SolidList ) )
        for j in m.SolidList:
            m.x[z,j,t].set_value( m.q[z,j,t].value/m.qT[z,t].value )
        if z != 0:
            for i in m.rxn_idx:
                calculate_variable_from_constraint(m.k[z,i,t],m.eq_r1[z,i,t])
            calculate_variable_from_constraint(m.X[z,t],m.eq_r2[z,t])
            m.X_term[z,t].set_value( (1-m.X[z,t].value)**(2/3) )
            calculate_variable_from_constraint(m.X_term[z,t],m.eq_r3[z,t])
            for i in m.rxn_idx:
                calculate_variable_from_constraint(m.r_gen[z,i,t],m.eq_r4[z,i,t])



    

    #with open('post_alg_update.txt','w') as f:
    #    m.display(ostream=f)

def integrate(fs,t): 
    # function to integrate the model one time step, from t
    m = fs.MB_fuel
    with open('dPdz.txt','w') as f:
        m.dPdz_disc_eq.pprint(ostream=f)
        m.dPdz_disc_eq.display(ostream=f)

    alg_update(fs,t)






def main():
    # main routine
    print('executing main routine')

if __name__ == '__main__':
    main()
