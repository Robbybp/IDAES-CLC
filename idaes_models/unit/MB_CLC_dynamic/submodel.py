# - initialize to steady state
# - apply perturbation to inputs
# - fix all algebraic variables (requires preciese knowledge of the 
#   algebraic variables) to their steady state values
# - deactivate all algebraic equations
# - solve single finite element (should be a dynamic problem, because inputs reach 
#   differential variables through boundary conditions)
# - algebraic update (?) (based on algebraic equations, which are deactivated) 
# - (re-fix algebraic variables... does setting value unfix variables?)
# - unfix/activate "first wave" of algebraic variables/equations
#   (could be as little as a single constraint/variable)
#   (may or may not change trajectory of differential variables)


# what happens when you have a problem, solved at a point, then
# add a constraint/variable (add a dimension), and try to solve for that 
# varaible?
# (what properties of the constraint set, projected into the dimension of the new 
# variable, guarantee convergence ??? (if any... - convexity, likely) )

# implicit euler step: algebraic variable at t0 is free (alg_update)
# algebraic variable at t1 is the only degree of freedom that has been added,
# but its initialization is incorrect, and when it gets recalculated to 
# satisfy its algebraic equation, the differential equations are no longer 
# satisfied 

# why do I expect this problem to be any easier to solve than just solving the
# the problem all at once ?????
# ^ no clear reason, yet, but it's at least something to try
from __future__ import division
from __future__ import print_function
from pyomo.environ import *
from pyomo.dae import *
from pyomo.core.base.constraint import SimpleConstraint
from pyomo.core.base.var import SimpleVar, IndexedVar
import pdb

def find_algebraic_variables(fs,diff_vars,time_derivatives,inputs,disturbances,geometry):
    # input: lists containing the (global) names of the differential variables,
    #        input variables, and disturbance variables (and geometric parameters, which are for some reason 
    #        treated as variables)
    # output: list of the algebraic variables (or algebraic variable names?)

    m = fs.MB_fuel
    alg_vars = []
    for var in m.component_objects(Var):
        if var.name not in diff_vars and var.name not in inputs and \
        var.name not in disturbances and var.name not in geometry and var.name not in time_derivatives:
            alg_vars.append(var)

    return alg_vars

def get_alg_var_data(alg_var):
    alg_var_data = []
    for var in alg_var:
        for index in var:
            if index == None:
                if not var.fixed:
                    alg_var_data.append(var)
                break
            else:
                if not var[index].fixed:
                    alg_var_data.append(var[index])

    return alg_var_data

def make_alg_var_const_map(fs):
    m = fs.MB_fuel 
    a = {}
    # vs
    #pdb.set_trace()
    #if var.parent_component() == m.vs:
    #    return m.eq_f6[var.index()]
    #for index in m.vs: a[ m.vs[index].name ] = m.eq_f6[index]
    #
    #for t in m.t:
    ## Boundary condition constraints:
    ## ... all of these are overwritten
    #    a[ m.MW_vap[0,t].name ] = m.eq_p5[0,t]
    #    a[ m.P[0,t].name ] = m.eq_e1[0,t]
    #    a[ m.rho_vap[0,t].name ] = m.eq_p6[0,t]
    #    a[ m.vg_in[t].name ] = m.eq_f5[t] 
    #    a[ m.vg[0,t].name ] = m.eq_e2[0,t]

    #    for j in m.GasList:
    #        a[ m.G_flux[0,j,t].name ] = m.eq_f1[j,t]
    #        a[ m.F[0,j,t].name ] = m.eq_c1[0,j,t] 
    #        a[ m.Cg[0,j,t].name ] = m.eq_c4[0,j,t] 
    #    for j in m.SolidList:
    #        a[ m.S_flux[1,j,t].name ] = m.eq_f2[j,t]
    #        a[ m.q[1,j,t].name ] = m.eq_c13[1,j,t] 

    # Gas composition and relevant properties:
    for index in m.CgT: a[ m.CgT[index].name ] = m.eq_c18[index]
    for index in m.eq_c6: a[ m.y[index].name ] = m.eq_c6[index]
    # ^ constraint indices used because y(z0) is specified by a recently added constraint
    for index in m.ytot: a[ m.ytot[index].name ] = m.eq_c7[index]
    for index in m.mu_vap: a[ m.mu_vap[index].name ] = m.eq_p8[index]
    for index in m.MW_vap: a[ m.MW_vap[index].name ] = m.eq_p5[index]
    for index in m.P: a[ m.P[index].name ] = m.eq_e1[index] 
    for index in m.rho_vap: a[ m.rho_vap[index].name ] = m.eq_p6[index]
    for index in m.dPdz_disc_eq: a[ m.dPdz[index].name ] = m.dPdz_disc_eq[index]
    for index in m.vg: a[ m.vg[index].name ] = m.eq_e2[index]
    for index in m.F: a[ m.F[index].name ] = m.eq_c4[index]
    for index in m.Ftotal: a[ m.Ftotal[index].name ] = m.eq_c3[index]
    for index in m.G_flux: a[ m.G_flux[index].name ] = m.eq_c1[index]
    for index in m.Gas_M: a[ m.Gas_M[index].name ] = m.eq_c2[index]
    for index in m.umf: a[ m.umf[index].name ] = m.eq_e3[index]
    for index in m.v_diff: a[ m.v_diff[index].name ] = m.eq_e4[index]
    
    # Solid composition and relevant properties
    for index in m.qT: a[ m.qT[index].name ] = m.eq_c14[index]
    for index in m.S_flux: a[ m.S_flux[index].name ] = m.eq_c13[index]
    for index in m.Solid_M: a[ m.Solid_M[index].name ] = m.eq_c9[index]
    for index in m.Solid_M_total: a[ m.Solid_M_total[index].name ] = m.eq_c10[index]
    for index in m.Solid_F: a[ m.Solid_F[index].name ] = m.eq_c11[index]
    for index in m.Solid_F_total: a[ m.Solid_F_total[index].name ] = m.eq_c12[index]

    # Reaction variables:
    for index in m.x: a[ m.x[index].name ] = m.eq_c15[index]
    for index in m.xtot: a[ m.xtot[index].name] = m.eq_c16[index]
    for index in m.X: a[ m.X[index].name ] = m.eq_r2[index]
    for index in m.eq_r1: a[ m.k[index].name ] = m.eq_r1[index]
    for index in m.eq_r3: a[ m.X_term[index].name ] = m.eq_r3[index]
    for index in m.eq_r4: a[ m.r_gen[index].name ] = m.eq_r4[index]
    for index in m.eq_r5: a[ m.rs[index].name ] = m.eq_r5[index]
    for index in m.eq_r6: a[ m.rg[index].name ] = m.eq_r6[index]
    for index in m.eq_b3: a[ m.Ctrans[index].name ] = m.eq_b3[index]
    for index in m.eq_b4: a[ m.qtrans[index].name ] = m.eq_b4[index]

    for index in m.X_gas: a[ m.X_gas[index].name ] = m.eq_c8[index]
    for index in m.X_OC: a[m.X_OC[index].name ] = m.eq_c17[index]

    # Thermal properties:
    for index in m.k_vap: a[ m.k_vap[index].name ] = m.eq_p16[index]
    for index in m.cp_vap: a[ m.cp_vap[index].name ] = m.eq_p10[index]
    for index in m.cp_gas: a[ m.cp_gas[index].name ] = m.eq_p11[index]
    for index in m.Pr: a[ m.Pr[index].name ] = m.eq_g3[index]
    for index in m.Rep: a[ m.Rep[index].name ] = m.eq_g1[index]
    for index in m.Nu: a[ m.Nu[index].name ] = m.eq_g6[index]
    for index in m.hf: a[ m.hf[index].name ] = m.eq_g9[index]
    for index in m.Tg_GS: a[ m.Tg_GS[index].name ] = m.eq_d3[index]
    for index in m.eq_p2: a[ m.DH_rxn_s[index].name ] = m.eq_p2[index]
    for index in m.eq_d7: a[ m.Ts_dHr[index].name ] = m.eq_d7[index]
    for index in m.cp_sol: a[ m.cp_sol[index].name ] = m.eq_p4[index]
    for index in m.cv_vap: a[ m.cv_vap[index].name ] = m.eq_p12[index]
    for index in m.k_cpcv: a[ m.k_cpcv[index].name ] = m.eq_p13[index]

    # Energy flux:
    for index in m.Gh_flux: a[ m.Gh_flux[index].name ] = m.eq_d2[index]
    for index in m.Sh_flux: a[ m.Sh_flux[index].name ] = m.eq_d11[index]

    # z-derivatives:
    for index in m.dG_fluxdz_disc_eq: a[ m.dG_fluxdz[index].name ] = m.dG_fluxdz_disc_eq[index]
    for index in m.dS_fluxdz_disc_eq: a[ m.dS_fluxdz[index].name ] = m.dS_fluxdz_disc_eq[index]
    for index in m.dGh_fluxdz_disc_eq: a[ m.dGh_fluxdz[index].name ] = m.dGh_fluxdz_disc_eq[index]
    for index in m.dSh_fluxdz_disc_eq: a[ m.dSh_fluxdz[index].name ] = m.dSh_fluxdz_disc_eq[index]

    # Outlet conditions:
    for index in m.Gas_Out_P: a[ m.Gas_Out_P[index].name ] = m.eq_f7[index]
    for index in m.Solid_Out_M: a[ m.Solid_Out_M[index].name ] = m.eq_f8[index]
    for index in m.Solid_Out_Ts: a[ m.Solid_Out_Ts[index].name ] = m.eq_f9[index]
    for index in m.Solid_Out_x: a[ m.Solid_Out_x[index].name ] = m.eq_f10[index]


    return a

def fix_z0(fs):
    # fix, at z=0, the variables of the model whose specifying constraints do not exist here
    # these are reaction-based algebraic constraints and z-derivatives
    m = fs.MB_fuel
    for t in m.t:
        for i in m.rxn_idx:
            m.k[0,i,t].fix(0)
            m.r_gen[0,i,t].fix(0)
        m.X_term[0,t].fix(0)
        m.Ts_dHr[0,t].fix(0)
        m.DH_rxn_s[0,t].fix(0)
        m.dGh_fluxdz[0,t].fix(0)
        m.dSh_fluxdz[0,t].fix(0)
        m.dPdz[0,t].fix(0)
        m.dldz[0].fix(m.L)
        # ^ unclear where this line should go
        for j in m.SolidList:
            m.rs[0,j,t].fix(0)
            m.qtrans[0,j,t].fix(0)
            m.dS_fluxdz[0,j,t].fix(0)
        for j in m.GasList:
            m.rg[0,j,t].fix(0)
            m.Ctrans[0,j,t].fix(0)
            m.dG_fluxdz[0,j,t].fix(0)
        
def unfix_z0(fs):
    m = fs.MB_fuel
    for t in m.t:
        for i in m.rxn_idx:
            m.k[0,i,t].unfix()
            m.r_gen[0,i,t].unfix()
        m.X_term[0,t].unfix()
        m.Ts_dHr[0,t].unfix()
        m.DH_rxn_s[0,t].unfix()
        m.dGh_fluxdz[0,t].unfix()
        m.dSh_fluxdz[0,t].unfix()
        m.dPdz[0,t].unfix()
        m.dldz[0].unfix()
        for j in m.SolidList:
            m.rs[0,j,t].unfix()
            m.qtrans[0,j,t].unfix()
            m.dS_fluxdz[0,j,t].unfix()
        for j in m.GasList:
            m.rg[0,j,t].unfix()
            m.Ctrans[0,j,t].unfix()
            m.dG_fluxdz[0,j,t].unfix()
