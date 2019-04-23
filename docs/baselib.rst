IDAES Standard Unit Models Library
==================================

This document provides a guide to the contents of the IDAES Standard Unit Model Library and how to use the models contained therein. The purpose of this library is to provide a set of simple models for common unit operations, which can either be used directly in process flowsheets, or as the basis for more detailed models of these unit operations.

All models within this library are designed to work within the standard IDAES framework and make use of Property classes for calculating the necessary physical, transport and reaction properties.

A list of the standard unit models available within the library is provided below. All standard models can be found in (idaes_models.core.unit.standard).

.. contents:: Contents 
    :depth: 3

Mixers and Splitters
--------------------
This section describes unit models with the purpose of combining or splitting material flows.

Mixer
^^^^^
*TBA*

Splitter
^^^^^^^^
This is a basic model for dividing a material stream into an arbitrary number of streams.

.. module:: idaes_models.unit.standard.splitter

.. autoclass:: _Splitter
    :special-members: __init__
    :members:

This class creates a model for a splitter unit with **Nout** outlets. The inlet flow is divided between each outlet based on a set of split fractions, **sfrac[Nout]**, with all outlets having the same state as the inlet stream.

**Assumptions**

1. All streams have the same state (T,P,y)

**Property Requirements**

The splitter model requires no property calculations, however a property package is still required to provide a component list for the unit.

**Degrees of Freedom**

Total degrees of freedom are *3 + number of components + (Nout-1)*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), along with all but one of the split fractions (**sfrac[]**) (the final split fraction is determined by the requirement that sum(sfrac) = 1). Alternatively, it is possible to specify the flowrate, pressure or composition at any outlet instead of the inlet.

**Initialization**

The default initialization routine allows the user to specify a state (**F**, **T**, **P** and **y**) at which they wish to initialize a splitter unit (specified by **blk**). This input is ignored for any state which is specified as *fixed* in the inlet stream. If the user does not provide a state and the inlet stream is not fixed, then the initialization defaults to *F=1.0*, *T=298.15*, *P=101325* and *y=1/no. comps*.

The initialization routine solves for the outlet flowrates based on the provided split fractions, at the specified temperature, pressure and composition.

Separator
^^^^^^^^^
*TBA*

Heaters and Coolers
-------------------
This section describes unit models with the purpose of implementing temperature changes.

Heater/Cooler
^^^^^^^^^^^^^
This is a basic model for a heat transfer unit with a defined heat duty or temperature differential.

.. module:: idaes_models.unit.standard.heater

.. autoclass:: _Heater
    :special-members: __init__
    :members:

The model utilizes the property package specified in **prop_lib** (with optional arguments **prop_lib_args**) to calculate the specific enthalpy of the inlet and outlet streams as a function of temperature, pressure and composition. An energy balance is then used to determine the heat duty of the unit.

**Assumptions**

1. Outlet composition equal to inlet composition
2. No pressure drop across unit

**Property Requirements**

The heater model requires only the total specific enthalpy (**h_mix**) at both the inlet and the outlet as a function of state (T,P,y).

**Degrees of Freedom**

Total degrees of freedom are *4 + number of components*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), along with one piece of information regarding the heat duty of the system (either heat duty (**Q**), temperature differential (**deltaT**) or outlet temperature (**Out_T**)). Alternatively, it is possible to specify the flowrate, pressure or composition at the outlet instead of the inlet.

**Initialization**

The default initialization routine allows the user to specify a state (**F**, **T**, **P** and **y**) at which they wish to initialize a heater unit (specified by **blk**). This input is ignored for any state which is specified as *fixed* in the inlet stream. If the user does not provide a state and the inlet stream is not fixed, then the initialization defaults to *F=1.0*, *T=298.15*, *P=101325* and *y=1/no. comps*.

The initialization routine first initializes the property blocks for the inlet and outlet, and then solves the mass and energy balances for the unit.

Heat Exchanger
^^^^^^^^^^^^^^
This is a basic model for a heat exchanger unit with a defined heat transfer coefficient and area.

.. module:: idaes_models.unit.standard.heat_exchanger

.. autoclass:: _HX
    :special-members: __init__
    :members:

This model consists of two Heater units (user may specify models to use with **heater_1** and **heater_2**) coupled via an energy balance with the constraint :math:`Q = U A dT`. Direction of flow (e.g. co-current or counter-current) can be specified using **flow** and the method used to calculate *dT* using **dT_method** (e.g. log-mean temperature difference). Properties for each stream are calculated using the property packages specified in **prop_lib_1** and **prop_lib_2** (with optional arguments **prop_lib_1_args** and **prop_lib_2_args**). If only **prop_lib_1** is specified then **prop_lib_2** defaults to the same package (along with any arguments give).

Each side of the heat exchanger is a separate Pyomo Block, named **Side_1** and **Side_2** respectively.

Equations are written such that either side of the heat exchanger can be the hot or cold side.

**Assumptions**

Assumptions in the heat exchanger unit model are the same as those in the selected Heater models.

**Property Requirements**

Property requires in the heat exchanger unit model are the same as those in the selected Heater models.

**Degrees of Freedom**

Degrees of freedom are *6 + number of components side 1 + number of components side 2 + 2* (plus any additional degrees of freedom from heater models). These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions) of both inlet streams, plus two variables selected from the heat transfer coefficient (**U**), heat transfer area (**A**) and outlet temperatures (**Side_1.Out_T** and **Side_2.Out_T**).

Alternatively, it is possible to specify the flowrate, pressure or composition at the outlets instead of the inlets.

**Initialization**

The default initialization routine allows the user to specify a state (**F**, **T**, **P** and **y**) for each side of a heat exchanger (specified by **blk**) at which they wish to initialize the heater sub-blocks. This input is ignored for any state which is specified as *fixed* in the inlet stream. If the user does not provide a state and the inlet stream is not fixed, then the initialization defaults to *F=1.0*, *T=298.15*, *P=101325* and *y=1/no. comps*.

The initialization routine first calls the initialization routine for each heater block at the specified initial guesses. After this, the cross-over temperature is estimated based on the specific enthalpies of each stream and used to get initial values for the heat duty on each side of the heat exchanger. Finally, the two heater units are coupled and the temperature driving force constraint enforced.


Pressure Changers
-----------------
This section describes unit models with the purpose of changing the pressure of a material flow.

Simple Pressure Changer
^^^^^^^^^^^^^^^^^^^^^^^
This is a basic model for a simple pressure changer unit.

.. module:: idaes_models.unit.standard.pressure_changer

.. autoclass:: _Simple
    :special-members: __init__
    :members:

This model consists of a basic pressure changer unit operation, where work (**W**) is equal to the enthalpy difference between the inlet and outlet streams with a specified outlet temperature. This model is intended primarily as a basis for more advanced models.

**Assumptions**

1. Outlet temperature known
2. :math:`W = F*(h_2-h_1)`

**Property Requirements**

The simple pressure changer model requires only the total specific enthalpy (**h_mix**) at both the inlet and the outlet as a function of state (T,P,y).

**Degrees of Freedom**

Total degrees of freedom are *5 + number of components*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), along with the outlet temperature (**Out_T**) and either the outlet pressure (**Out_P**) or pressure ratio (**ratioP**). Alternatively, it is possible to specify the flowrate, pressure or composition at the outlet instead of the inlet.

**Initialization**

The default initialization routine allows the user to specify a state (**F**, **T**, **P** and **y**) at which they wish to initialize a pressure changer unit (specified by **blk**). This input is ignored for any state which is specified as *fixed* in the inlet stream. If the user does not provide a state and the inlet stream is not fixed, then the initialization defaults to *F=1.0*, *T=298.15*, *P=101325* and *y=1/no. comps*.

The initialization routine first initializes the property blocks for the inlet and outlet, and then solves the mass and energy balances for the unit.

Isentropic Pressure Changer
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a model for an isentropic pressure changer unit.

.. autoclass:: _Isentropic
    :special-members: __init__
    :members:

This model builds upon the simple pressure changer model by adding constraints for isentropic operation and mechanical efficiency. Model can be built to represent a compressor/pump with **compressor** = *True*, or as an expander with **compressor** = *False*. This flag changes whether mechanical work (**W**) is greater than (compressor) or less than (expander) the isentropic work (**W_isen**). The mechanical efficiency (**eta**) of the unit is assumed to be fixed.

**Assumptions**

1. Isentropic operation
2. Fixed mechanical efficiency

**Property Requirements**

The isentropic pressure changer model requires a property package that provides the total specific enthalpy and total specific entropy (**h_mix** and **s_mix**) as a function of state (T,P,y).

**Degrees of Freedom**

Total degrees of freedom are *5 + number of components*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), along with the pressure ratio (**ratioP**) and mechanical efficiency (**eta**). 

Alternatively, it is possible to specify the flowrate, pressure or composition at the outlet instead of the inlet, as well as work in place of pressure ratio.

**Initialization**

The default initialization routine allows the user to specify a state (**F**, **T**, **P** and **y**) at which they wish to initialize an isentropic pressure changer unit (specified by **blk**). This input is ignored for any state which is specified as *fixed* in the inlet stream. If the user does not provide a state and the inlet stream is not fixed, then the initialization defaults to *F=1.0*, *T=298.15*, *P=101325* and *y=1/no. comps*.

The initialization routine first initializes the property blocks for the inlet and outlet, and then solves the mass and energy balances for the unit.

Compressor
^^^^^^^^^^
This is a model for an isentropic compressor unit.

.. module:: idaes_models.unit.standard.compressor

.. autoclass:: _Compressor
    :special-members: __init__
    :members:

The compressor model class is a wrapper for the isentropic pressure changer model, which sets the **compressor** flag to *True*. See the documentation for the isentropic pressure changer model (above) for more details.

Expander
^^^^^^^^
This is a model for an isentropic expander unit.

.. module:: idaes_models.unit.standard.expander

.. autoclass:: _Expander
    :special-members: __init__
    :members:

The expander model class is a wrapper for the isentropic pressure changer model, which sets the **compressor** flag to *False*. See the documentation for the isentropic pressure changer model (above) for more details.

Phase Equilibrium Models
------------------------
This section describes unit models with the purpose of performing phase equilibrium operations.

Flash Unit
^^^^^^^^^^

This is a model for basic two-phase flash unit.

.. module:: idaes_models.unit.standard.flash

.. autoclass:: _Flash
    :special-members: __init__
    :members:

This model utilizes property packages capable of vapor-liquid equilibrium calculations to perform a single stage, vapor-liquid flash operation. The actual phase equilibrium calculations are performed by the property package, whilst the unit model handles only the overall mass and energy balances.

This model can be used for all types of flash operations (e.g. TP, TH, PH etc).

**Assumptions**

1. System is at equilibrium

**Property Requirements**

The flash unit operation requires a property package capable of performing phase equilibrium calculations, and makes use of the following properties:

* Vapor fraction - **vf**
* Liquid phase composition - **y_liq**
* Vapor phase composition - **y_vap**
* Total Specific Enthalpy - **h_mix**

**Degrees of Freedom**

Total degrees of freedom are *5 + number of components*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), along with two variables chosen from the unit heat duty (**Q**) and the state of the outlet (such as temperature (**Out_T**) and pressure (**Out_P**)). 

Alternatively, it is possible to specify the flowrate, pressure or composition at the outlet instead of the inlet.

**Initialization**

The default initialization routine allows the user to specify a state (**F**, **T**, **P** and **y**) at which they wish to initialize a basic reactor unit (specified by **blk**). This input is ignored for any state which is specified as *fixed* in the inlet stream. If the user does not provide a state and the inlet stream is not fixed, then the initialization defaults to *F=1.0*, *T=298.15*, *P=101325* and *y=1/no. comps*.

The initialization routine first initializes the property blocks for the inlet and outlet, and then solves the mass and energy balances for the unit.

Tray Distillation Column
^^^^^^^^^^^^^^^^^^^^^^^^
*TBA*

Reactor Models
--------------
This section describes unit models for common reaction operations.

Stoichiometric Reactor
^^^^^^^^^^^^^^^^^^^^^^
This is a model for stoichiometric or yield reactor unit.

.. module:: idaes_models.unit.standard.rstoic

.. autoclass:: _RStoic
    :special-members: __init__
    :members:

This is a model for a basic reaction operation where the user specifies either the extents of reaction (**x_rxn**) or composition of the outlet stream (**Out_y**) (i.e. as either a stoichiometric or yield reactor). 

**Assumptions**

1. Outlet temperature and pressure known

**Property Requirements**

The stoichiometric reactor unit model requires a property package which provides:

* a set of valid reactions (**rxn_idx**),
* an array of stoichiometric coefficients for all possible reactions (**stoic**, indexed by reaction and components),
* total specific enthalpy (**h_mix**) as a function of state (T,P,y).

**Degrees of Freedom**

Total degrees of freedom are *5 + number of components + number of reactions*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), the outlet pressure (**Out_P**), either the outlet temperature (**Out_T**) or heat duty (**Q**) and either one extent of reaction (**x_rxn**) or outlet composition (**y_out**) for possible each reaction. 

**Initialization**

The default initialization routine allows the user to specify a state (**F**, **T**, **P** and **y**) at which they wish to initialize a pressure changer unit (specified by **blk**). This input is ignored for any state which is specified as *fixed* in the inlet stream. If the user does not provide a state and the inlet stream is not fixed, then the initialization defaults to *T=298.15*, *P=101325* and *y=1/no. comps*.

The initialization routine first initializes the property blocks for the inlet and outlet, and then solves the mass and energy balances for the unit.

CSTR
^^^^
This is a model for well mixed tank reactor (CSTR) unit.

.. module:: idaes_models.unit.standard.rcstr

.. autoclass:: _RCSTR
    :special-members: __init__
    :members:

This is a model for a well mixed reactor (CSTR), which builds upon the stoichiometric reactor unit model. The CSTR model adds constraints for the extents of reaction based upon the volume of the reactor vessel (**V**) and the rate of reaction calculated by a reaction property package.

**Assumptions**

1. Perfectly mixed tank
2. Temperature and pressure of reaction vessel known

**Property Requirements**

The CTR unit model requires a property package which provides:

* a set of valid reactions (**rxn_idx**),
* an array of stoichiometric coefficients for all possible reactions (**stoic**, indexed by reaction and components),
* rates of reaction (**rate_rxn**) for all reactions as a function of state,
* total specific enthalpy (**h_mix**) as a function of state (T,P,y).

**Degrees of Freedom**

Total degrees of freedom are *6 + number of components*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), plus the reactor volume (**V**), pressure (**Out_P**) and either the reactor temperature (**Out_P**) or heat duty (**Q**). 

**Initialization**

The default initialization routine is the same as for the stoichiometric reactor unit model (see above).

Equilibrium Reactor
^^^^^^^^^^^^^^^^^^^
This is a model for equilibrium reactor unit using equilibrium constants.

.. module:: idaes_models.unit.standard.requil

.. autoclass:: _REquil
    :special-members: __init__
    :members:

This is a model for an equilibrium reactor, which builds upon the stoichiometric reactor unit model. The equilibrium model adds constraints to solve for a state at which all rates of reaction (provided by a reaction property package) are equal to zero.

**Assumptions**

1. Chemical equilibrium - all rates of reaction are zero
2. Temperature and pressure of reaction vessel known

**Property Requirements**

The equilibrium reactor unit model requires a property package which provides:

* a set of valid reactions (**rxn_idx**),
* an array of stoichiometric coefficients for all possible reactions (**stoic**, indexed by reaction and components),
* rates of reactions (**rate_rxn**) for all reactions as a function of state,
* total specific enthalpy (**h_mix**) as a function of state (T,P,y).

**Degrees of Freedom**

Total degrees of freedom are *5 + number of components*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), the outlet pressure (**Out_P**) and either the outlet temperature (**Out_T**) or heat duty (**Q**).

**Initialization**

The default initialization routine is the same as for the stoichiometric reactor unit model (see above).

Gibbs Reactor
^^^^^^^^^^^^^
This is a model for equilibrium reactor unit using Gibbs energy minimization.

.. module:: idaes_models.unit.standard.rgibbs

.. autoclass:: _RGibbs
    :special-members: __init__
    :members:

This is a model for an equilibrium reactor based on minimization of Gibbs energy. This is achieved using Lagrange multipliers and the pure component Gibbs energies for all components in the system. Log flowrates are used to improve scaling of small flowrates.

**Assumptions**

1. Chemical equilibrium - minimized Gibbs energy of system
2. Temperature and pressure of reaction vessel known

**Property Requirements**

The Gibbs reactor unit model requires a property package which provides:

* a list of all chemical elements in the system (**elem**),
* elemental composition of all chemical species present (**elem_comp**),
* pure component specific Gibbs energies for all chemical species (**g_pc**) as a function of state (T,P),
* total specific enthalpy (**h_mix**) as a function of state (T,P,y).

**Degrees of Freedom**

Total degrees of freedom are *5 + number of components*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), the outlet pressure (**Out_P**) and either the outlet temperature (**Out_T**) or heat duty (**Q**).

**Initialization**

The default initialization routine allows the user to specify a state (**F**, **T**, **P** and **y**) at which they wish to initialize a Gibbs reactor unit (specified by **blk**). This input is ignored for any state which is specified as *fixed* in the inlet stream. If the user does not provide a state and the inlet stream is not fixed, then the initialization defaults to *F=1.0*, *T=298.15*, *P=101325* and *y=1/no. comps*.

The initialization routine first initializes the property blocks for the inlet and outlet, and then solves the mass and energy balances for the unit.

PFR
^^^
This is a model for plug flow reactor (PFR) unit.

.. module:: idaes_models.unit.standard.rpfr

.. autoclass:: _RPFR
    :special-members: __init__
    :members:

This is a model for a plug flow reactor (PFR) unit. This model uses a normalized axial domain (0 to 1) and makes use of the Pyomo DAE toolbox for discretizing the domain and calculating derivatives.

This model allows for side inputs for both material (**F_side**, **y_side** and **Fh_side**) and heat (**Q**), as well as pressure drop along the reactor (**deltaP**). These are *fixed* to zero by default, but users who wish to make use of these features may *unfix* these.

**Assumptions**

1. Plug flow
2. No pressure drop along reactor (can be easily replaced)

**Property Requirements**

The PFR unit model requires a property package which provides:

* a set of valid reactions (**rxn_idx**),
* an array of stoichiometric coefficients for all possible reactions (**stoic**, indexed by reaction and components),
* rates of reaction (**rate_rxn**) for all reactions as a function of state,
* total specific enthalpy (**h_mix**) as a function of state (T,P,y).

**Degrees of Freedom**

Total degrees of freedom are *5 + number of components*. These would normally be the feed flowrate (**In_F**), temperature (**In_T**), pressure (**In_P**) and composition (**In_y**) (mole fractions), plus the reactor length (**Lb**), and either the reactor area (**A**) or volume (**V**). 

**Initialization**

The default initialization routine allows the user to specify a state (**F**, **T**, **P** and **y**) at which they wish to initialize a PFR unit (specified by **blk**). This input is ignored for any state which is specified as *fixed* in the inlet stream. If the user does not provide a state and the inlet stream is not fixed, then the initialization defaults to *F=1.0*, *T=298.15*, *P=101325* and *y=1/no. comps*.

The initialization routine first initializes the property blocks for the entire domain, and then sequentially solves the mass balances, energy balances and pressure drop (if present) for the unit.

