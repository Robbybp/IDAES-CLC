General Cubic Equation of State Model
=====================================

This page describes the General Cubic Equation of State model developed as part of the IDAES standard model library.

.. contents:: Contents 
    :depth: 3

Reference
---------

Once this is published, a reference to the paper should be given here.

Introduction
------------

The IDAES General Cubic Equation of State model is a general framework for solving any cubic equation of state of the form:

.. math::
    0 = Z^3 - (1+B-uB)Z^2 + (A+[w-u]B^2-uB)Z - AB - wB^2 - wB^3

This form includes many of the common cubic equations of state, such as the Peng-Robinson and Soave-Redlich-Kwong equations.

The model takes the state of a fluid (e.g. temperature (**T**), pressure (**P**) and mole fractions (**y_mix**)) and calculates thermodynamic and phase equilibrium properties for the fluid. The model supports two phase (vapor-liquid equilibrium) and uses departure functions to calculate the other state variables of the system (including both individual phase and total mixture properties).

The model also allow the user to specify that their system has only a single phase, which suppresses the phase equilibrium calculations and simplifying the system of equations.

**IMPORTANT NOTE**

In order to ensure model equations are not singular, the phase equilibrium calculations must be solved in the two phase region. For conditions outside the two-phase region, the general cubic equation of state model locates the nearest point within the two phase region, and solves the phase equilibrium at that point (see Details of Model for more information). Due to this, users should be aware that properties for non-existent phases are solved at a state different to the specified system state, and property values will reflect this.

Properties Available in Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following properties are currently supported by the general cubic equation of state model.

* Compressibility Factor - **Z_liq**, **Z_vap**
* Critical Pressure - **Pc_mix**
* Critical Temperature - **Tc_mix**
* Fugacity Coefficients - **phi_liq**, **phi_vap**
* Molar Volume - **V_liq**, **V_vap**, **V_mix**
* Molecular Weight - **mw_liq**, **mw_vap**, **mw_mix**
* Phase Equilibrium - **y_liq**, **y_vap**
* Reduced Pressure - **Pr**
* Reduced Temperature - **Tr**
* Specific Enthalpy - **h_liq**, **h_vap**, **h_mix**
* Specific Entropy - **s_liq**, **s_vap**, **s_mix**
* Specific Gibbs Energy - **g_liq**, **g_vap**, **g_mix**
* Vapor fraction - **vf**

The following cubic equations of state are supported by the general cubic equation of state model.

* Peng-Robinson - **PR**
* Soave-Redlich-Kwong - **SRK**

Structure of the Model
----------------------

The IDAES general cubic equation of state model consists of three parts:

1. The core Cubic_EoS.py module, which contains most of the constraints necessary to solve the phase equilibrium and thermodynamic properties,
2. A user provided parameter file, which contains the necessary pure component parameters and information needed to specify the system of interest, and
3. An external function library, which contains functions to find the roots of a cubic equation.

Details on each of these components is given below, along with  information on how to use them.

Cubic EoS Module
^^^^^^^^^^^^^^^^

The PropPack class in the Cubic_EoS.py file contains the core components of the general cubic equation of state model.

.. module:: idaes_models.properties.physical.Cubic_EoS

.. autoclass:: PropPack
    :special-members: __init__
    :members:

The PropPack class is not meant to be used as a stand-alone model - it is designed to be called by a user defined parameter class via inheritance (see below). When called, the PropPack class expects to find certain objects and methods defined by the user parameter class, which are used in constructing property package.

The Cubic_EoS class also allows the user to specify a number of arguments that affect how the resulting Pyomo Block is constructed. These are:

* **plib** specifies the location of the external function call library to use for finding the roots of the cubic equation. By default, the model looks in the same directory as the Cubic_EoS.py file.
* **phase** indicates whether the model should consider the full two phase system or instead assume only a single phase is present. Specifying a single phase (**V** or **L**) results in a Pyomo Block which does not contain the equations for phase equilibrium. This greatly simplifies the model equations, at the expense of ignoring potential phase equilibria.
* **mix_props** specifies whether the model should calculate properties for the total mixture in addition to phase specific properties. Most models in the IDAES Standard Unit Model Library will expect the model to provide mixture properties.

Cubic Root Finder Function
^^^^^^^^^^^^^^^^^^^^^^^^^^

**NEED MORE DOCUMENTATION HERE (or a link)**

The cubic root finder function is an Pyomo ExternalFunction call, which is used to solve for the roots of the cubic equation, and to provide the necessary partial derivatives to the solver. The advantage of this approach is that it allows for simple analytical methods to be used to find the roots of the equation of state, rather than relying on inequality constraints numerical methods to solve for the roots.

The IDAES model library contains a default root finding function, which can be found in (idaes_models.core.properties.physical.cubic_eos). This directory contains a c-file containing the code for the root finder and a Makefile for compiling the code on Linux systems. This library must first be compiled before it can be used, which can be done using the following instructions:

1. Copy the files **phys_prop.c** and **Makefile** to the desired install directory (we recommend placing these in the same directory as Cubic_EoS.py).
2. In a shell prompt, navigate to the install directory.
3. Execute the make command.

Assuming your system has a compatible compiler, this should compile the phys_prop.c file, creating a file named **phys_prop.so**. This file can then be used for for the general cubic equation of state by pointing the **plib** argument to the **phys_prop.so** file. If you compiled phys_prop.so in the same directory as Cubic_EoS.py, then it is not necessary to specify plib, as the model will search for phys_prop.so in the same directory by default.

Parameter Files
^^^^^^^^^^^^^^^

The parameter files provides the user with the means to specify the necessary parameters for defining their system of interest. The parameter file needs to contain the following elements in order to make use of the general cubic equation of state model:

1. A class named **_PropPack** which inherits from the *PropPack* class in the Cubic_EoS.py file, along with a **ProcBlock** decorator (see documentation elsewhere regarding ProcBlock decorator).
2. **comp** - a Python list containing names for each chemical species to be considered.
3. A method named **_thermo_params** which contains:
    * **EoS_ID** - an identifier for which equation of state to use (currently supports *PR* and *SRK*)
    * **elem** - A list of all chemical elements present in the system
    * **elem_comp** - a dictionary indexed by **comp**, specifying the elemental composition of each species
    * **mw_pc** - molecular weights for each component species
    * **Tc** - critical temperatures for each component species
    * **Pc** - critical pressures for each component species
    * **omega** - Pitzer acentricity factors for each component species
    * **kappa** - binary interaction parameters for each species pair (dependent on the chosen equation of state)
    * **Antoine** - coefficients for the Antoine equation (3 parameter version) for each species. **N.B.** units for the Antoine coefficients are in **bar** and **K**, rather than SI. This is the form Antoine coefficients are most commonly given in.

4. A method named **_ideal_pc_props** which contains:
    * A Pyomo constraint indexed by **comp** which calculates the pure component, ideal gas heat capacity (**cp_ig_pc**) for each species as a function of state (T,P)
    * A Pyomo constraint indexed by **comp** which calculates the pure component, ideal gas specific enthalpy (**h_ig_pc**) for each species as a function of state (T,P)
    * A Pyomo constraint indexed by **comp** which calculates the pure component, ideal gas specific entropy (**s_ig_pc**) for each species as a function of state (T,P)
    * Any parameters necessary for the above constraints

An example of a property class is given in (idaes_models.core.properties.physical.PR_ASU), and discussed below. Users are encouraged to use this example as a template for developing their own parameter classes.

**Example Property Class**

This is an example of a property class for using the Peng-Robinson equation of state to calculate the properties of mixtures of Ar, N2 and O2.

.. module:: idaes_models.properties.physical.PR_ASU

.. autoclass:: _PropPack
    :special-members: __init__
    :members:

As the _PropPack class in the PR_ASU.py file inherits from the PropPack class in Cubic_EoS.py, all the arguments are the same.

Using the Model
---------------

In order to make use of the general cubic equation of state model with their desired parameter file, users should import their parameter file into their unit models and call the **PropPack** class when creating property blocks in their models. For example:

``self.prop = ParameterFile.PropPack(name='Prop', parent=self, [arguments])``

**N.B.** The PropPack class is a metaclass created from the **_PropPack** class via the **ProcBlock** decorator.

**Degrees of Freedom**

The general cubic equation of state model has *2 + number of components* degrees of freedom (assuming all parameters are specified). These would normally be the composition of the entire fluid mixture (**y_mix**) and two other state variables (most often temperature (**T**) and pressure (**P**)).

**Model Parameters**

The general cubic equation of state has the following parameters that users may wish to adjust for their circumstances.

* **P_ref** - thermodynamic reference pressure (default *P_ref = 101,325 Pa*)
* **T_ref** - thermodynamic reference temperature (default *T_ref = 298.15 K*)
* **eps1**, **eps2** - smoothing factors for smooth maximum operations (default *eps1 = 1e-3* and *eps2 = 5e-4*)

**Initialization**

The default initialization routine allows the user to specify a state (**T**, **P** and **y_mix**) at which they wish to initialize a property block (specified by **blk**).  This input is ignored for any state which is specified as *fixed*. If the user does not provide a state and the variable is not fixed, then the initialization defaults to *T=298.15*, *P=101325* and *y_mix=1/no. comps*.

Optionally, users may also specify initial guesses for the liquid and vapor phase compositions (**y_liq** and **y_vap**) if they have them.

The initialization routine uses the following sequence:

1. Initialize pure component properties
2. Initialize intermediate parameters for phase equilibrium
3. Initialize complementarity variables for phase equilibrium
4. Solve phase equilibrium
5. Calculate liquid phase properties
6. Calculate vapor phase properties
7. Calculate mixture properties

Some of these steps may be skipped if the user chooses to limit the model to consider only certain phases.

Creating New Parameter Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users are encouraged to make use of the example parameter file when developing their own parameter files.

Details of Model
----------------

Developing a general form for solving phase equilibrium using a cubic equation of state model is complicated by three factors:

1. Identifying and solving for the correct roots of the cubic equation,
2. Dealing with conditions where one root is imaginary (single phase regions), and
3. Dealing with singularities that arise when there is only a single root (supercritical region).

A number of approaches have been developed by various groups to address theses issues, and these can be found in the open literature. For the IDAES General Cubic Equation of State Model, the approach taken to deal with these factors is outlines below.

Solving the Cubic Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^

A number of methods for identifying and solving the cubic equation within an equation oriented framework have been proposed by different groups, such as directly solving the cubic equation subject to inequality constraints to identify the correct root. However, there also exist a number of simple analytical methods for solving for the roots of cubic equations from which it is a simple matter to proceduraly select the desired root. It is also relatively easy to analytically calculate partial derivatives for the the roots of the cubic equation of state as a function of the A and B terms commonly contained therein.

By making use of Pyomo's link to the AMPL Solver Library (ASL) and its ability to make calls to external functions contained in compiled library file, the IDAES General Cubic Equation of State model makes use of these analytical and procedural approaches to solving for the correct roots of the cubic equation of state. This greatly simplifies the problem being passes to the solver (whilst retaining partial derivative information) and removes the need for inequality constraints. This also allows for "pseudo-roots" to be returned in cases where one of the desired roots has an imaginary component, which prevents infeasible solutions from arising if the solver moves outside the two-phase envelope.

Disappearing Phases
^^^^^^^^^^^^^^^^^^^

The next problem that needs to be addressed is how to deal with states where only one phase can exist. In these cases, there is no valid solution to the cubic equation for the non-existent phase, and thus it is not possible to compute properties for this phase either. This creates problems in a general equation oriented framework, as it is necessary to calculate values for all variable at all states.

To address this issue, the IDAES General Cubic Equation of State Model instead solves the phase equilibrium at a temperature **Te** which satisfies the following conditions:

1. If both phase exist at state temperature **T**, then :math:`Te = T`, otherwise
2. **Te** is equal to the nearest point on the boundary of the two-phase region at the current state pressure, **P**.

By solving the phase equilibrium problem at **Te** rather than **T**, this ensures that there will always be a valid solution for both phases. For states where **T** falls within the two-phase region this results in the expected solution. For situations where **T** lies outside the two-phase region, by using a **Te** that lies on the boundary of the two-phase region, we allow for a valid solution for both phases, with a negligible amount of the non-existent phase. **NOTE** that this means that the General Cubic Equation of State Model will always predict some amount of material existing in both phases, however one phase will often be negligible.

Supercritical Conditions
^^^^^^^^^^^^^^^^^^^^^^^^

In conditions at or above the critical point, the solution of the phase equilibrium problem is complicated by the fact that the roots for both phases are equal, and thus vapor and liquid phase are indistinguishable. If not addressed, this results in a singularity as the fugacities for each component in each phase are equal, and thus the amount of vapor and liquid cannot be determined.

In the IDAES General Cubic Equation of State Model, this has been addressed by adding to the complementarity constraints for the disappearing phase problem to enforce that the vapor fraction in the supercritical region will be 0 (i.e. all liquid).

Implementation
^^^^^^^^^^^^^^

To implement the above solutions, a new variable, **beta**, defined such that :math:`beta = T/Te`, and a system of complementarity constraints is applied to enforce the following logical constraints (**vf** is the vapor phase fraction and **lf** is the liquid phase fraction:

1. If **P > Pc_mix**, the system will always be liquid and **vf = 0**
2. If **P < Pc_mix** and **T > Tc_mix**, the system will always be vapor and **vf = 1**
3. Otherwise (if **P < Pc_mix** and **T < Tc_mix**), the system is potentially a two phase mixture and one of the following must be true:
    a. If **beta > 1** then the system is in the liquid only region and **lf = 0**
    b. If **beta < 1** then the system is in the vapor only region and **vf = 0**
    c. Otherwise **beta = 0**

These can be written using MAX operations to express the logical tests, which are implemented using smooth maximum approximations. The smoothing parameters **eps1** and **eps2** are used to applying the smoothing to the maximum operators.

Due to the use of smooth maximum operations, there is a degree of error introduced into the calculations which users need to be aware of. Most notably, this manifests in that the vapor and liquid fractions (**vf** and **lf**) will never be equal to zero. In some cases, the non-existent phase may even report a negligible negative value. Users need to be aware of this behavior, and include robustness checks in their models to deal with any issues that may arise because of these.

Phase Properties
^^^^^^^^^^^^^^^^

As phase properties are dependent on a valid solution to the cubic equation of state for the phase in question, it is also necessary to use **Te** as the state when calculating phase properties for the non-existent phase. Users need to be aware of this, as this affects the values that properties will take.

Assumptions
^^^^^^^^^^^

1. Supercritical fluid is represented as liquid phase.




