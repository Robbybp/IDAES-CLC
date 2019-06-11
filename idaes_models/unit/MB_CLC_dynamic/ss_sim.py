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

from pyomo.environ import value
from pyomo.opt import SolverFactory

import time

from idaes_models.core import FlowsheetModel, ProcBlock
import SS_MB_CLC as MB_CLC_fuel

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
        Args:
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
        """           
        # Create unit model for fuel reactor
        self.MB_fuel = MB_CLC_fuel.MB(
                parent=self,
                dae_method = 'OCLR',
                press_drop = 'Ergun',
                fe_set = fe_set,
                ncp = 3)

        
def setInputs(fs):    
    # ===== Fuel Reactor ===== 
    # Gas phase inlet conditions
    fs.MB_fuel.Gas_In_F.fix(128.20513)   # mol/s
    fs.MB_fuel.Gas_In_P.fix(2.00)   # bar *estimated. min pressure to overcome pressure drop
    fs.MB_fuel.Gas_In_Tg.fix(293.15)     # K
    fs.MB_fuel.Gas_In_y['CO2'].fix(0.02499)
    fs.MB_fuel.Gas_In_y['H2O'].fix(0.00001)
    fs.MB_fuel.Gas_In_y['CH4'].fix(0.975) 
       
    # Solid phase inlet conditions
    fs.MB_fuel.Solid_In_M.fix(591.4) #479.011) # kg/s
    fs.MB_fuel.Solid_In_Ts.fix(1183.15)      # K
    fs.MB_fuel.Solid_In_x['Fe2O3'].fix(0.44999)
    fs.MB_fuel.Solid_In_x['Fe3O4'].fix(1e-5)
    fs.MB_fuel.Solid_In_x['Al2O3'].fix(0.55)
    
    # Bed characteristics
    fs.MB_fuel.Dr.fix(6.5) # m
    fs.MB_fuel.L.fix(5) # m
    fs.MB_fuel.eps.fix(0.4) # (-)   
    
    
def print_summary_fuel_reactor(fs):
    """
    Print some key results.
    """
    print("\nResults:")
    print("==========================================")
    print("---Moving Bed Fuel Reactor---")       
    
    print("\nInlet gas: ", 
              "\nCO2: ", value(fs.MB_fuel.F[0,'CO2']), "mol/s",
              "\nH20: ", value(fs.MB_fuel.F[0,'H2O']), "mol/s",
              "\nCH4: ", value(fs.MB_fuel.F[0,'CH4']), "mol/s",
              "\nCO2: ", value(fs.MB_fuel.Gas_M[0,'CO2']), "kg/s",
              "\nH20: ", value(fs.MB_fuel.Gas_M[0,'H2O']), "kg/s",
              "\nCH4: ", value(fs.MB_fuel.Gas_M[0,'CH4']), "kg/s")
    print("\nOutlet gas: ", 
              "\nCO2: ", value(fs.MB_fuel.F[1,'CO2']), "mol/s",
              "\nH20: ", value(fs.MB_fuel.F[1,'H2O']), "mol/s", 
              "\nCH4: ", value(fs.MB_fuel.F[1,'CH4']), "mol/s",
              "\nCO2: ", value(fs.MB_fuel.Gas_M[1,'CO2']), "kg/s",
              "\nH20: ", value(fs.MB_fuel.Gas_M[1,'H2O']), "kg/s", 
              "\nCH4: ", value(fs.MB_fuel.Gas_M[1,'CH4']), "kg/s")
    print("\nInlet solids: ", 
              "\nFe2O3: ", value(fs.MB_fuel.Solid_F[1,'Fe2O3']), "mol/s",
              "\nFe3O4: ", value(fs.MB_fuel.Solid_F[1,'Fe3O4']), "mol/s", 
              "\nAl: ", value(fs.MB_fuel.Solid_F[1,'Al2O3']), "mol/s",
              "\nFe2O3: ", value(fs.MB_fuel.Solid_M[1,'Fe2O3']), "kg/s",
              "\nFe3O4: ", value(fs.MB_fuel.Solid_M[1,'Fe3O4']), "kg/s", 
              "\nAl: ", value(fs.MB_fuel.Solid_M[1,'Al2O3']), "kg/s")
    print("\nOutlet solids: ", 
              "\nFe2O3: ", value(fs.MB_fuel.Solid_F[0,'Fe2O3']), "mol/s",
              "\nFe3O4: ", value(fs.MB_fuel.Solid_F[0,'Fe3O4']), "mol/s", 
              "\nAl: ", value(fs.MB_fuel.Solid_F[0,'Al2O3']), "mol/s",
              "\nFe2O3: ", value(fs.MB_fuel.Solid_M[0,'Fe2O3']), "kg/s",
              "\nFe3O4: ", value(fs.MB_fuel.Solid_M[0,'Fe3O4']), "kg/s", 
              "\nAl: ", value(fs.MB_fuel.Solid_M[0,'Al2O3']), "kg/s") 
    
    print("\nGas inlet velocity: ", value(fs.MB_fuel.vg[0]), "m/s")
    print("Gas outlet velocity: ", value(fs.MB_fuel.vg[1]), "m/s")
    print("Solids velocity: ", value(fs.MB_fuel.vs), "m/s")    
    
    print("\nHeat of reaction @ z=0: ", 
              value(fs.MB_fuel.DH_rxn_s[0]), "J/(mol reaction)")
    print("Heat of reaction @ z=1: ", 
              value(fs.MB_fuel.DH_rxn_s[1]), "J/(mol reaction)")
    
    print("\nCH4 conversion: ", value(fs.MB_fuel.X_gas)*100, " %")
    print("Fe2O3 conversion: ", value(fs.MB_fuel.X_OC)*100, " %")
    
    print('\nPressure @inlet: ', value(fs.MB_fuel.P[0]))
    print('Pressure @outlet: ', value(fs.MB_fuel.Gas_Out_P))
    
    print("\nReactor bed height:", value(fs.MB_fuel.L), " m")
    print("Reactor bed diameter:", value(fs.MB_fuel.Dr), " m")
#    print("Refractory wall thickness", value(fs.MB.refractory_th), " m")
    
    print("\nInlet gas flow:", value(fs.MB_fuel.Gas_In_F), " mol/s")
    print("Outlet gas flow:", value(fs.MB_fuel.Ftotal[1]), " mol/s")
    print("Inlet solids flow:", value(fs.MB_fuel.Solid_In_M), " kg/s")
    print("Outlet solids flow:", value(fs.MB_fuel.Solid_Out_M), " kg/s")
    print("Inlet solids temperature:", value(fs.MB_fuel.Solid_In_Ts), " K")
    print("Outlet solids temperature:", value(fs.MB_fuel.Solid_Out_Ts), " K")
    
    print("Inlet gas temperature:", value(fs.MB_fuel.Tg[0]), " K")
    print("Outlet gas temperature:", value(fs.MB_fuel.Tg[1]), " K")    
    
    print("\nInlet solid mass fractions: ", 
              "\nFe2O3: ", value(fs.MB_fuel.x[1,'Fe2O3']),
              "\nFe3O4: ", value(fs.MB_fuel.x[1,'Fe3O4']), 
              "\nAl2O3: ", value(fs.MB_fuel.x[1,'Al2O3']))
    print("Outlet solid mass fractions: ", 
              "\nFe2O3: ", value(fs.MB_fuel.x[0,'Fe2O3']),
              "\nFe3O4: ", value(fs.MB_fuel.x[0,'Fe3O4']), 
              "\nAl2O3: ", value(fs.MB_fuel.x[0,'Al2O3'])) 
    
def results_plot_fuel_reactor(self):
    """
    Plot some key results.
    """
    
    import matplotlib.pyplot as plt    

    # Total pressure profile
    P = []
    for z in self.MB_fuel.z:
        P.append(value(self.MB_fuel.P[z]))
    fig_P = plt.figure(1)
    plt.plot(self.MB_fuel.z, P)
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Total Pressure [bar]")       

    # Temperature profile
    Tg = []
    Ts = []
#    Tw = []
    for z in self.MB_fuel.z:
        Tg.append(value(self.MB_fuel.Tg[z] - 273.15))
        Ts.append(value(self.MB_fuel.Ts[z] - 273.15))
#        Tw.append(value(self.MB_fuel.Tw[z]))
    fig_T = plt.figure(2)
    plt.plot(self.MB_fuel.z, Tg, label='Tg')
    plt.plot(self.MB_fuel.z, Ts, label='Ts')
#    plt.plot(self.MB_fuel.z, Tw, label='Tw')
    plt.legend(loc=0,ncol=2)
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Temperature [C]") 
    
    # Superficial gas velocity and minimum fluidization velocity
    vg = []
    umf = []
    for z in self.MB_fuel.z:
        vg.append(value(self.MB_fuel.vg[z]))
        umf.append(value(self.MB_fuel.umf[z]))
    fig_vg = plt.figure(3)
    plt.plot(self.MB_fuel.z, vg, label='vg')
    plt.plot(self.MB_fuel.z, umf, label='umf')
    plt.legend(loc=0,ncol=2)
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Superficial gas velocity [m/s]")
    
   # Gas components molar flow rate
    for j in self.MB_fuel.GasList:
        F = []
        for z in self.MB_fuel.z:
            F.append(value(self.MB_fuel.F[z,j]))
        fig_F = plt.figure(4)
        plt.plot(self.MB_fuel.z, F, label=j)
    plt.legend(loc=0,ncol=len(self.MB_fuel.GasList))
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Gas component molar flow rate, F [mol/s]")  
    
    # Bulk gas phase total molar flow rate
    Ftotal = []
    for z in self.MB_fuel.z:
        Ftotal.append(value(self.MB_fuel.Ftotal[z]))
    fig_Ftotal = plt.figure(5)
    plt.plot(self.MB_fuel.z, Ftotal)
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Total molar gas flow rate [mol/s]")  

    # Solid components mass flow rate
    for j in self.MB_fuel.SolidList:
        M = []
        for z in self.MB_fuel.z:
            M.append(value(self.MB_fuel.Solid_M[z,j]))
        fig_M = plt.figure(6)
        plt.plot(self.MB_fuel.z, M, label=j)
    plt.legend(loc=0,ncol=len(self.MB_fuel.SolidList))
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Solid components mass flow rate [kg/s]")
    
     # Bulk solid phase total molar flow rate
    Mtotal = []
    for z in self.MB_fuel.z:
        Mtotal.append(value(self.MB_fuel.Solid_M_total[z]))
    fig_Mtotal = plt.figure(7)
    plt.plot(self.MB_fuel.z, Mtotal)
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Solid total mass flow rate [kg/s]")        
    
    # Gas phase concentrations
    for j in self.MB_fuel.GasList:
        Cg = []
        for z in self.MB_fuel.z:
            Cg.append(value(self.MB_fuel.Cg[z,j]))
        fig_Cg = plt.figure(8)
        plt.plot(self.MB_fuel.z, Cg, label=j)
    plt.legend(loc=0,ncol=len(self.MB_fuel.GasList))
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Concentration [mol/m3]")       
    
    # Gas phase mole fractions
    for j in self.MB_fuel.GasList:
        y = []
        for z in self.MB_fuel.z:
            y.append(value(self.MB_fuel.y[z,j]))
        fig_y = plt.figure(9)
        plt.plot(self.MB_fuel.z, y, label=j)
    plt.legend(loc=0,ncol=len(self.MB_fuel.GasList))
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("y [-]")  
    
    # Solid phase mass fractions
    for j in self.MB_fuel.SolidList:
        x = []
        for z in self.MB_fuel.z:
            x.append(value(self.MB_fuel.x[z,j]))
        fig_x = plt.figure(10)
        plt.plot(self.MB_fuel.z, x, label=j)
    plt.legend(loc=0,ncol=len(self.MB_fuel.SolidList))
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("x [-]")  

    # Total mass fraction
    xtot = []
    for z in self.MB_fuel.z:
        xtot.append(value(self.MB_fuel.xtot[z]))
    fig_xtot = plt.figure(11)
    plt.plot(self.MB_fuel.z, xtot)
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Total mass fraction [-]") 
    
    # # Gas mix density
    # rhog = []
    # for z in self.MB_fuel.z:
    #     rhog.append(value(self.MB_fuel.rho_vap[z]))
    # fig_rhog = plt.figure(23)
    # plt.plot(self.MB_fuel.z, rhog)
    # plt.grid()
    # plt.xlabel("Bed height [-]")
    # plt.ylabel("Gas mix density [kg/m3]") 
               
    # Fe conversion
    X_Fe = []
    for z in self.MB_fuel.z:
        X_Fe.append(value(self.MB_fuel.X[z])*100)
    fig_X_Fe = plt.figure(13)
    plt.plot(self.MB_fuel.z, X_Fe)
    plt.grid()
    plt.xlabel("Bed height [-]")
    plt.ylabel("Fraction of metal oxide converted [%]") 
                       
       
def main():
    """
    Make the flowsheet object and solve
    """
    flowsheet = Flowsheet(name='MB_Model') 
    
    # Fix variables
    setInputs(flowsheet) 

    ts = time.time() 

    mb = flowsheet.MB_fuel
    
    # Initialize fuel reactor
    flowsheet.MB_fuel._initialize(outlvl=1,
                              optarg={"tol"            : 1e-8,
                                      "max_cpu_time"   : 600,
                                      "print_level"    : 1,
                                      "halt_on_ampl_error": 'yes'})        
       
    # Create a solver
    opt = SolverFactory('ipopt')
    opt.options = {'tol': 1e-8,
                   'linear_solver'  : 'ma27',
                   'bound_push': 1e-8,
                   'max_cpu_time': 600,
                   'print_level': 5}
    
    results = opt.solve(flowsheet,tee=True,symbolic_solver_labels=False,
                            keepfiles=False)

    #flowsheet.MB_fuel.Solid_In_M.fix(691.4)
    #flowsheet.MB_fuel.Gas_In_y['CO2'].fix(0.03999)
    #flowsheet.MB_fuel.Gas_In_y['H2O'].fix(0.00001)
    #flowsheet.MB_fuel.Gas_In_y['CH4'].fix(0.96)



    #results = opt.solve(flowsheet,tee=True,symbolic_solver_labels=False,
    #                        keepfiles=False)
    
    
    print("\n")
    print("----------------------------------------------------------")
    print('Total simulation time: ', value(time.time() - ts), " s")
    print("----------------------------------------------------------")

    
    # Print some variables 
    #print_summary_fuel_reactor(flowsheet) 

    # Plot some variables 
    #results_plot_fuel_reactor(flowsheet) 

    with open('ss_flowsheet.txt','w') as f:
        flowsheet.display(ostream=f)


    dt_Gflux_CO2    = []
    dt_Gflux_H2O    = []
    dt_Gflux_CH4    = []
    dt_Sflux_Fe2O3  = []
    dt_Sflux_Fe3O4  = []
    dt_Sflux_Al2O3  = []
    dt_Ctrans_CO2   = []
    dt_Ctrans_H2O   = []
    dt_Ctrans_CH4   = []
    dt_qtrans_Fe2O3 = []
    dt_qtrans_Fe3O4 = []
    dt_qtrans_Al2O3 = []
    dt_Ghflux       = []
    dt_Ts           = []
    dt_TgGS         = []
    dt_TsGS         = []
    dt_vg           = []
    dt_vs           = []

#    for z in mb.z.get_finite_elements():
#        if z != mb.z.first() and z != mb.z.last():
#
#            dt_Gflux_CO2.append( (mb.Cg[z,'CO2'].value-mb.Cg[prev,'CO2'].value)/\
#                    (mb.G_flux[z,'CO2'].value-mb.G_flux[prev,'CO2'].value) \
#                    *(z-prev)*mb.eps.value*mb.L.value /(z-prev))
#
#            dt_Gflux_H2O.append( (mb.Cg[z,'H2O'].value-mb.Cg[prev,'H2O'].value)/\
#                    (mb.G_flux[z,'H2O'].value-mb.G_flux[prev,'H2O'].value) \
#                    *(z-prev)*mb.eps.value*mb.L.value /(z-prev))
#
#            dt_Gflux_CH4.append( (mb.Cg[z,'CH4'].value-mb.Cg[prev,'CH4'].value)/\
#                    (mb.G_flux[z,'CH4'].value-mb.G_flux[prev,'CH4'].value) \
#                    *(z-prev)*mb.eps.value*mb.L.value /(z-prev))
#
#            dt_Ctrans_CO2.append( (mb.Cg[z,'CO2'].value-mb.Cg[prev,'CO2'].value)/\
#                    (mb.Ctrans[z,'CO2'].value)* \
#                    #-mv.Ctrans[prev,'CO2'].value)*\
#                    mb.eps.value/(1-mb.eps.value) /(z-prev))
#
#            dt_Ctrans_H2O.append( (mb.Cg[z,'H2O'].value-mb.Cg[prev,'H2O'].value)/\
#                    (mb.Ctrans[z,'H2O'].value)* \
#                    #-mv.Ctrans[prev,'H2O'].value)*\
#                    mb.eps.value/(1-mb.eps.value) /(z-prev))
#
#            dt_Ctrans_CH4.append( (mb.Cg[z,'CH4'].value-mb.Cg[prev,'CH4'].value)/\
#                    (mb.Ctrans[z,'CH4'].value)* \
#                    #-mv.Ctrans[prev,'CH4'].value)*\
#                    mb.eps.value/(1-mb.eps.value) /(z-prev))
#
#            dt_Sflux_Fe2O3.append( (mb.q[z,'Fe2O3'].value-mb.q[prev,'Fe2O3'].value)/\
#                    (mb.S_flux[z,'Fe2O3'].value-mb.S_flux[prev,'Fe2O3'].value)*\
#                    (z-prev)/(1-mb.eps.value)*mb.L.value /(z-prev))
#
#            dt_Sflux_Fe3O4.append( (mb.q[z,'Fe3O4'].value-mb.q[prev,'Fe3O4'].value)/\
#                    (mb.S_flux[z,'Fe3O4'].value-mb.S_flux[prev,'Fe3O4'].value)*\
#                    (z-prev)/(1-mb.eps.value)*mb.L.value /(z-prev))
#
#            dt_Sflux_Al2O3.append( (mb.q[z,'Al2O3'].value-mb.q[prev,'Al2O3'].value)/\
#                    (mb.S_flux[z,'Al2O3'].value-mb.S_flux[prev,'Al2O3'].value)*\
#                    (z-prev)/(1-mb.eps.value)*mb.L.value /(z-prev))
#
#            dt_qtrans_Fe2O3.append( (mb.q[z,'Fe2O3'].value-mb.q[prev,'Fe2O3'].value)/\
#                    (mb.qtrans[z,'Fe2O3'].value )/(z-prev)) 
#                    #-mb.qtrans[prev,'Fe2O3'].value) )
#
#            dt_qtrans_Fe3O4.append( (mb.q[z,'Fe3O4'].value-mb.q[prev,'Fe3O4'].value)/\
#                    (mb.qtrans[z,'Fe3O4'].value )/(z-prev)) 
#                    #-mb.qtrans[prev,'Fe3O4'].value) )
#
#            dt_qtrans_Al2O3.append( (mb.q[z,'Fe3O4'].value-mb.q[prev,'Fe3O4'].value)/\
#                    (mb.qtrans[z,'Fe3O4'].value )/(z-prev)) 
#                    #-mb.qtrans[prev,'Fe3O4'].value) )
#
#            dt_Ghflux.append( (mb.Tg[z].value-mb.Tg[prev].value)/\
#                (mb.Gh_flux[z].value-mb.Gh_flux[prev].value)* (z-prev)* mb.eps.value*\
#                mb.L.value* mb.rho_vap[z].value* mb.cp_gas[z].value /(z-prev)) 
#
#            dt_Ts.append( (z-prev)*(1-mb.eps.value)*mb.L.value/mb.vs.value /(z-prev))
#
#            dt_TgGS.append( (mb.Tg[z].value - mb.Tg[prev].value)/\
#                    mb.Tg_GS[z].value* mb.eps.value* mb.rho_vap[z].value* mb.cp_gas[z].value 
#                    /(z-prev))
#    
#            dt_TsGS.append( (mb.Ts[z].value - mb.Ts[prev].value)/\
#                    mb.Tg_GS[z].value* (1-mb.eps.value)* mb.rho_sol.value* mb.cp_sol[z].value*1e-3 
#                    /(z-prev))
#    
#            dt_vg.append( mb.L.value*(z-prev)/mb.vg[z].value /(z-prev))
#    
#            dt_vs.append( mb.L.value*(z-prev)/mb.vs.value /(z-prev))
#
#        prev = z
#
#    with open('dt.txt','w') as f:
#        f.write('dt_Gflux_CO2\t')
#        for t in dt_Gflux_CO2:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_Gflux_H2O\t')
#        for t in dt_Gflux_H2O:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_Gflux_CH4\t')    
#        for t in dt_Gflux_CH4:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_Sflux_Fe2O3\t') 
#        for t in dt_Sflux_Fe2O3:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_Sflux_Fe3O4\t')  
#        for t in dt_Sflux_Fe3O4:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_Sflux_Al2O3\t')  
#        for t in dt_Sflux_Al2O3:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_Ctrans_CO2\t')  
#        for t in dt_Ctrans_CO2:
#            f.write('%1.3f'%t +'\t')
#            
#        f.write('\ndt_Ctrans_H2O\t')   
#        for t in dt_Ctrans_H2O:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_Ctrans_CH4\t')   
#        for t in dt_Ctrans_CH4:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_qtrans_Fe2O3\t') 
#        for t in dt_qtrans_Fe2O3:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_qtrans_Fe3O4\t') 
#        for t in dt_qtrans_Fe3O4:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_qtrans_Al2O3\t') 
#        for t in dt_qtrans_Al2O3:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_Ghflux\t')      
#        for t in dt_Ghflux:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_Ts\t\t')           
#        for t in dt_Ts:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_TgGS\t\t')         
#        for t in dt_TgGS:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_TsGS\t\t')         
#        for t in dt_TsGS:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_vg\t\t')           
#        for t in dt_vg:
#            f.write('%1.3f'%t +'\t')
#
#        f.write('\ndt_vs\t\t')           
#        for t in dt_vs:
#            f.write('%1.3f'%t +'\t')


    # Store the flowsheet    
    return flowsheet
    
if __name__ == "__main__":
    flowsheet = main()         
