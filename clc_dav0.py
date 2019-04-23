"""
Flowsheet model for the BFB_CLC fuel reactor

Updated on Tue Aug 30 2018

@authors: cokoli           
"""
from __future__ import division
from __future__ import print_function
#import sys as sys

__author__ = "Chinedu Okoli"
__version__ = "1.0.0"

from pyomo.environ import value
from pyomo.opt import SolverFactory
import matplotlib.pyplot as plt

from idaes_models.core import FlowsheetModel, ProcBlock

from BFB_CLC_A import BFB_CLC # If BFB_CLC_A is not in the current working directory, a full path will be needed

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
        # Create set of finite elements for FR
        nfe = 10
        fe_a = 1/4.0
        fe_b = 0.2
        fe_set = [0, 0.004] # Custom grid size used for 1st element because of rapid rxn at gas entry
        for i in range(1,nfe+1):
            if i < nfe*fe_a:
                fe_set.append(i*fe_b/(nfe*fe_a))
            elif i == nfe:
                fe_set.append(1.0)
            else:
                fe_set.append(fe_b + (i-nfe*fe_a)*(1-fe_b)/(nfe*(1-fe_a)))

        # Create unit models for FR
        self.BFB_FR = BFB_CLC(parent=self,dae_method = 'OCLR',
                                fe_set = fe_set,
                                ncp = 3,
                                s_inlet = "Bottom",
                                s_outlet = "Overflow",
                                hx_type = "None",
                                prop_lib = 'FR_props',
                                hx_prop_lib = [],
                                vb_method = 'Davidson')

def setInputs(fs):
    #Fix some variables for the FR
    fs.BFB_FR.Dt.fix(6.5) # m
    fs.BFB_FR.Lb.fix(5) # m
    fs.BFB_FR.nor.fix(2500) # (-)   

    fs.BFB_FR.Gas_In_F.fix(272.81)   # mol/s 272.81 mol/s is the value based from Peltola et al
    fs.BFB_FR.Gas_In_P.fix(156000)   # Pa *estimated. min pressure to overcome pressure drop
    fs.BFB_FR.Gas_In_T.fix(550)        # K
    fs.BFB_FR.Gas_In_y['CO2'].fix(0.4772)
    fs.BFB_FR.Gas_In_y['H2O'].fix(0.0646)
    fs.BFB_FR.Gas_In_y['CH4'].fix(0.4582)
    
    fs.BFB_FR.Solid_In_F.fix(1422) # kg/s
    fs.BFB_FR.Solid_In_T.fix(1186)      # K
    fs.BFB_FR.Solid_In_x['Fe2O3'].fix(0.45)
    fs.BFB_FR.Solid_In_x['Fe3O4'].fix(0.00)
    fs.BFB_FR.Solid_In_x['Al2O3'].fix(0.55)


def print_summary(fs):  
    print('FR_Inputs')
    fs.BFB_FR.Gas_In_F.display()
    fs.BFB_FR.Gas_In_T.display()
    fs.BFB_FR.Gas_In_P.display()
    fs.BFB_FR.Gas_In_y.display()
    fs.BFB_FR.Solid_In_F.display()
    fs.BFB_FR.Solid_In_T.display()
    fs.BFB_FR.Solid_In_x.display()
    
    print ()
    print('FR_Outputs')
    fs.BFB_FR.Gas_Out_F.display()
    fs.BFB_FR.Gas_Out_T.display()
    fs.BFB_FR.Gas_Out_P.display()
    fs.BFB_FR.Gas_Out_y.display()
    fs.BFB_FR.Solid_Out_F.display()
    fs.BFB_FR.Solid_Out_T.display()
    fs.BFB_FR.Solid_Out_x.display()
 
    print()
    print('FR_Size')
    fs.BFB_FR.Lb.display()
    fs.BFB_FR.Dt.display() 
  
#    # Calculations for material and energy balance tolerance
##    removal = {}
    mbal_tol = {}
    OC_conv = {} 
               
    OC_conv['Fe2O3'] = (fs.BFB_FR.Solid_In_F.value \
                                * fs.BFB_FR.Solid_In_x['Fe2O3'].value \
                            - fs.BFB_FR.Solid_Out_F.value \
                                * fs.BFB_FR.Solid_Out_x['Fe2O3'].value) \
                            / (fs.BFB_FR.Solid_In_F.value \
                                * fs.BFB_FR.Solid_In_x['Fe2O3'].value)
                         
    for j in fs.BFB_FR.GasList:
        mbal_tol[j] = (fs.BFB_FR.Gas_In_F.value
                            * fs.BFB_FR.Gas_In_y[j].value
                        - fs.BFB_FR.Gas_Out_F.value
                            * fs.BFB_FR.Gas_Out_y[j].value) \
                        / (fs.BFB_FR.Gas_In_F.value
                            * fs.BFB_FR.Gas_In_y[j].value)        

    ebal_gas_FR = value(fs.BFB_FR.Gas_In_F) \
                    * value(fs.BFB_FR.prop_b[0].h_vap) \
                - value(fs.BFB_FR.Gas_Out_F) \
                    * value(fs.BFB_FR.gas_prop_out.h_vap) 

    if fs.BFB_FR.s_outlet == 'Overflow':
        # Overflow conditions
        ebal_sol_FR = value(fs.BFB_FR.Solid_In_F)*value(fs.BFB_FR.sol_prop_f.h_sol) \
                    - value(fs.BFB_FR.Solid_Out_F) \
                        * value(fs.BFB_FR.prop_e[1].h_sol)
    else:
        # Underflow conditions
        ebal_sol_FR = value(fs.BFB_FR.Solid_In_F)*value(fs.BFB_FR.sol_prop_f.h_sol) \
                    - value(fs.BFB_FR.Solid_Out_F) \
                        * value(fs.BFB_FR.prop_e[0].h_sol)
    ebal_tol_FR = ebal_gas_FR + ebal_sol_FR + value(fs.BFB_FR.Q)

##    print('Removal:',removal)
    print()
    print('OC Conversion_FR:', OC_conv)
    print('Mass Balance Tolerance_FR:',mbal_tol)
    print()
    print('Energy Balance Tolerance_FR:',ebal_tol_FR)
    print('Energy balance gas_FR:',ebal_gas_FR)
    print('Energy balance solids_FR:', ebal_sol_FR)


#%% Output options

def results_plot_FR(self):
#    print('Fuel Reactor plots')
    Tge = []
    Tgc = []
    Tgb = []
    Tse = []
    Tsc = []
    Ge = []
    Gb = []
    cbt = []
    cct = []
    cet = []

    for i in self.BFB_FR.l_n:
        Tge.append(value(self.BFB_FR.Tge[i]))
        Tgc.append(value(self.BFB_FR.Tgc[i]))
        Tgb.append(value(self.BFB_FR.Tgb[i]))
        Tse.append(value(self.BFB_FR.Tse[i]))
        Tsc.append(value(self.BFB_FR.Tsc[i]))

        Ge.append(value(self.BFB_FR.Ge[i]))
        Gb.append(value(self.BFB_FR.Gb[i]))

        cbt.append(value(self.BFB_FR.cbt[i]))
        cct.append(value(self.BFB_FR.cct[i]))
        cet.append(value(self.BFB_FR.cet[i]))
        
    #Tray temperature profile
    plt.figure(1)
    plt.plot(self.BFB_FR.l_n, Tge, label='Tge')
    plt.plot(self.BFB_FR.l_n, Tgc, label='Tgc')
    plt.plot(self.BFB_FR.l_n, Tgb, label='Tgb')

    plt.legend(loc=9,ncol=2)
    plt.grid()
    plt.xlabel("Bed height")
    plt.ylabel("Gas temperatures in bed regions (K)")
    
    plt.figure(2)
    plt.plot(self.BFB_FR.l_n, Tse, label='Tse')
    plt.plot(self.BFB_FR.l_n, Tsc, label='Tsc')

    plt.legend(loc=9,ncol=3)
    plt.grid()
    plt.xlabel("Bed height")
    plt.ylabel("Solid temperatures in bed regions (K)")

    plt.figure(3)
    plt.plot(self.BFB_FR.l_n, Ge, label='Ge')
    plt.plot(self.BFB_FR.l_n, Gb, label='Gb')

    plt.legend(loc=9,ncol=2)
    plt.grid()
    plt.xlabel("Bed height")
    plt.ylabel("Gas flow (mol/s)")
    
    plt.figure(4)
    plt.plot(self.BFB_FR.l_n, cbt, label='cbt')
    plt.plot(self.BFB_FR.l_n, cct, label='cct')
    plt.plot(self.BFB_FR.l_n, cet, label='cet')

    plt.legend(loc=9,ncol=3)
    plt.grid()
    plt.xlabel("Bed height")
    plt.ylabel("gas conc. mol/s")
    
    #vapour phase mole composition
    for i in self.BFB_FR.GasList:
        y=[]
        for j in self.BFB_FR.l_n:
            y.append(value(self.BFB_FR.yb[i,j]))
        plt.figure(5)
        plt.plot(self.BFB_FR.l_n,y,label=i)
    plt.legend(loc=9,ncol=len(self.BFB_FR.GasList))
    plt.grid()
    plt.xlabel("Bed height")
    plt.ylabel("Gas bubble mole frac. (-)")
    
    #solid phase mass composition
    for i in self.BFB_FR.SolidList:
        x=[]
        for j in self.BFB_FR.l_n:
            x.append(value(self.BFB_FR.xe[i,j]))
        plt.figure(6)
        plt.plot(self.BFB_FR.l_n,x,label=i)
    plt.legend(loc=9,ncol=len(self.BFB_FR.SolidList))
    plt.grid()
    plt.xlabel("Bed height")
    plt.ylabel("Solid emulsion mass frac. (-)")
    
#%%
## Script to initialize, run the model and plot some results
if __name__ == "__main__":
    flowsheet = Flowsheet()
    
    setInputs(flowsheet)
    
    # Initialize models   
    print() 
    print("Initialize FR")
    
    flowsheet.BFB_FR._initialize(outlvl=0,
                              optarg={"tol"            : 1e-7,
                                      "max_cpu_time"   : 300,
                                      "print_level"    : 3,
                                      "bound_push" : 1e-8,
                                      "halt_on_ampl_error": 'yes'})
    print() 
    print("Simulation solve of FR")
    sopts = {"tol"            : 1e-7,
             "max_cpu_time"   : 300,
             "print_level"    : 3,
             "mu_init"        : 1e-3,
             "bound_push"     : 1e-5}
    opt = SolverFactory('ipopt') 
    results = opt.solve(flowsheet.BFB_FR,tee=True,keepfiles=False,options=sopts)
    

    print_summary(flowsheet)
    results_plot_FR(flowsheet)

# -----------------------------------------------------------------------------
