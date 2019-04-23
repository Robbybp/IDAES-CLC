IDAES Modeling Standards
========================

This is a document of the current standards for IDAES models as of 8th May 2017. This is still a work in progress, and is likely to change on a regular basis.

.. contents:: Contents 
    :depth: 3

Units and Reference State
-------------------------

The standard system of units for IDAES models is SI without prefixes (i.e. Pa, not kPa). The default thermodynamic reference state is 298.15 K and 101325 Pa.

For supercritical fluids, the standard approach is to consider the fluid as being part of the liquid phase (as it will be handled via pumps rather than compressors).

Model Classes
-------------

The current approach to developing models as part of the IDAES framework is to use callable classes to construct Pyomo Blocks, which can be connected together into larger flowsheets. These classes may be imported into a model and then called as part of a rule for a Pyomo Block object.

**Add documentation of ProcBlock decorator, etc., or a reference to such**

The current approach to naming these classes is:

* Unit Models - name of unit operation, or abbreviation thereof
* Property Packages - **PropPack**

Defining Component Lists
------------------------

*This section is very much a work in progress, based on my experiences whilst developing both the standard unit models, cubic EoS and BFB models (ALee). Feedback welcome*

All property packages should contain a list of components, **comp**, considered in the property package. This list of components should be created in such a way that it is accessible to unit models before the property package is constructed. One way to achieve this is to initially create the component list as a Python list (with the name **comp**) as part of the property package class, which can then be used to initialize a Pyomo Set. This allows any model using the property package to look up the **comp** Python list before constructing the property package class itself.

**Component Lists for Multiphase Systems**

For property packages involving multiphase systems, it is recommended that the property package contain component lists for each phase, as well as a full list of components in all phases (i.e. the union of all phase component lists). In these cases, the recommended naming convention is:

* Full List (Union) - **comp**
* Gas or Vapor Phase - **comp_vap**
* Liquid Phase - **comp_liq**
* Solid Phase - **comp_sol**

In situation where there are multiple phases of one type (such as two liquid phases), an integer should be added to the phase identifier, e.g. **_liq1**, **_liq2**. Integers should start at 1 and should be continuous.

When defining variables indexed by component list, modelers should default to indexing by the full (union) list unless it is obvious that a given property only exists for certain phases.

*End WIP*

Ports and Connecting Models
---------------------------

The standard definitions for material ports in IDAES are:

* Fluid ports - flowrate (**F**), temperature (**T**), pressure (**P**), mole fractions (**y**) and vapor fraction (**vf**)
* Solid ports - flowrate (**F**), temperature (**T**), pressure (**P**), mass fractions (**x**) and vapor fraction (**vf**)

In addition, energy ports are also supported, which contain a single energy flow (**E**).

Model developers should assume that all input ports are fully and consistently defined (e.g. that mole fractions sum to 1), so there should be no need to add constraints on the inlet flows. Similarly, developers should ensure that all outlet ports are fully and consistently defined to ensure downstream models do not receive inconsistent data.

Model Formatting
----------------

The section describes the recommended formatting for all models developed as part of the IDAES framework.

Headers and Meta-data
^^^^^^^^^^^^^^^^^^^^^

All models should include a header section which contains information on the model including:

* Model name,
* Model version number (see below),
* Model execution date (for simulation results),
* Model publication date,
* Model author,
* Pyomo version number,
* Modeling standards version number (see beginning of document).

* Any necessary licensing and disclaimer information (see below).
* Any additional information the modeler feels should be included.

Version Numbering
^^^^^^^^^^^^^^^^^

*TBD*

Licensing Information and Disclaimers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*TBD*

Coding Standard
^^^^^^^^^^^^^^^

All models and property packages produced as part of IDAES should conform to the PEP-8 standard.

Model Organization
^^^^^^^^^^^^^^^^^^

In addition to any constraints imposed by the IDAES framework, Python and Pyomo environments, models should be constructed in a logical fashion to aid other users in understanding the model. Model constraints should be grouped with similar constraints, and each grouping of constraints should be clearly commented. For unit operation models, the following grouping and ordering of constraints is recommended.

* Overall mass balances.
* Component mass balances (grouped by phase/region as necessary).
* Energy balances for each phase/region.
* Summation equations (sum of mole fractions, sum of component flows, etc.)
* Pressure relationships.
* Reaction and equilibrium constraints (if not part of a sub-model).
* Supporting equations (grouped with similar constraints as necessary).
* Any initialization constraints should be placed at the end of the model and clearly identified as such.

For property packages, it is recommended that all the equations necessary for calculating a given property be grouped together, clearly separated and identified by using comments.

Additionally, model developers are encouraged to consider breaking their model up into a number of smaller methods where this makes sense. This can facilitate modification of the code by allowing future users to inherit from the base model and selectively overloading sub-methods where desired.

Naming Conventions
^^^^^^^^^^^^^^^^^^

Currently, the only standards for naming conventions within IDAES models is a list of standard variable names to use for state variables and properties (see next section). For other elements of their models, modelers are encouraged to use a consistent approach to naming of elements (most notably variable and constraint names).

Many of the models currently in the IDAES model libraries identify constraints with the prefix **eq_**.

Commenting
^^^^^^^^^^

To help other modelers and users understand the how a model works, model builders are strongly encouraged to comment their code. It is suggested that every constraint should be commented with a description of the purpose of the constraint, and if possible/necessary a reference to a source or more detailed explanation. Any deviations from standard units or formatting should be clearly identified here. Any initialization procedures, or other procedures required to get the model to converge should be clearly commented and explained where they appear in the code. Additionally, modelers are strongly encouraged to add additional comments explaining how their model works to aid others in understanding the model.

Standard Variable Names
-----------------------

In order to maximize interoperability of models and code, a set of recommended names has been developed for common variables has been developed. Modelers are strongly recommended to make use of these standard names, as this will allow their models to interact with models from other developers using the standard with a minimum of overhead.

Standard variables names are primarily associated with state and property variables, and other variables which are commonly used to interact between units.

Phase Identification
^^^^^^^^^^^^^^^^^^^^

In many systems of interest, there will potentially be multiple phases present, each of which can have its own set of states and properties. To handle this, it is recommended that a phase identifier be appended to all property, and if necessary state, variables to unambiguously identify what each variable refers to. For variables referring to the average properties of a mixture (such as the lumped enthalpy of a two-phase stream), it is recommended that the variable be appended with the **_mix** identifier.

* Total Mixture - **_mix**
* Solids        - **_sol**
* Liquids       - **_liq**
* Vapor/Gases   - **_vap**

**Multiple Phases**

In situation where there are multiple phases of one type (such as two liquid phases), an integer should be added to the phase identifier, e.g. **_liq1**, **_liq2**. Integers should start at 1 and should be continuous.

**Pure Components**

In some situations it may also be necessary to specify properties for pure components within a system (a notable example is the need for pure component Gibbs energies for Gibbs energy minimization reactors). In these cases, the modifier **_pc** should be appended after the phase identifier (e.g. *h_liq_pc*).

State Variables
^^^^^^^^^^^^^^^

The standard variable names for defining the state of a stream are listed below. Mass fractions are recommended for solid phases, whilst mole fractions are recommended for fluid phases.

* Flowrate       - **F** [:math:`mol/s`]
* Temperature    - **T** [:math:`K`]
* Pressure       - **P** [:math:`Pa`]
* Mass Fractions - **x** [:math:`kg/kg`]
* Mole Fractions - **y** [:math:`mol/mol`]
* Vapor Fraction - **vf** [:math:`-`]

Thermodynamic Properties
^^^^^^^^^^^^^^^^^^^^^^^^

The standard variable names for thermodynamic properties are listed below.

* Critical Pressure - **Pc** [:math:`Pa`]
* Critical Temperature - **Tc** [:math:`K`]
* Density - **rho** [:math:`kg/m^3`]
* Fugacity - **f** [:math:`Pa`]
* Fugacity Coefficient - **phi** [:math:`-`]
* Henry's Constant - **henry** [*TBD*]
* Molecular Weight - **mw** [:math:`kg/mol`]
* Reduced Pressure - **Pr** [:math:`Pa`]
* Reduced Temperature - **Tr** [:math:`K`]
* Specific Heat Capacity (const. P) - **cp** [:math:`J/mol.K`]
* Specific Heat Capacity (const. V) - **cv** [:math:`J/mol.K`]
* Specific Heat of Reaction - **dH** [:math:`J/mol`]
* Specific Enthalpy - **h** [:math:`J/mol`]
* Specific Entropy  - **s** [:math:`J/mol.K`]
* Specific Gibbs Energy - **g** [:math:`J/mol`]
* Specific Helmholtz Energy - **a** [:math:`J/mol`]
* Specific Volume - **V** [:math:`m^3/mol`]

Transport Properties
^^^^^^^^^^^^^^^^^^^^

The standard variable names for transport properties are listed below.

* Diffusivity - **D** [:math:`m^2/s`]
* Surface Tension - **sigma** [:math:`N/m`]
* Thermal Conductivity - **kc** [:math:`W/m.K`]
* Viscosity (dynamic) - **mu** [:math:`kg/m.s` (:math:`Pa.s`)]
* Viscosity (kinematic) - **nu** [:math:`m^2/s`]

Reaction Properties
^^^^^^^^^^^^^^^^^^^

The standard variable names for properties associated with reaction rates and kinetics are listed below. Where possible, heats of reaction (or heats of formation) should be incorporated into the calculations for specific enthalpy to allow simplify the energy balances in unit models using property packages.

* List of identifiers for potential reactions - **rxn_idx**

* Activity - **act** [:math:`-`]
* Activity Coefficients - **gamma** [:math:`-`]
* Heat of Reaction - **dH_rxn** [:math:`J/mol`]
* Rate Constant - **k_rxn** [*varies*]
* Equilibrium Constant - **Keq** [:math:`-`]
* Rate of Reaction - **rate_rxn** [:math:`mol/m^3.s`]
* Stoichiometric coefficients - **stoic**

**Heterogeneous and Multiphase Systems**

For heterogeneous reaction systems, reaction properties should have an associated phase identifier to make it clear which phase is being used as the basis for the calculation. For example, if reaction rates are based on the solid phase (and thus the volume of solids), the **_sol** identifier should be appended where appropriate.

Solid Properties
^^^^^^^^^^^^^^^^

The standard variable names for properties associated with solids and particles are listed below.

* Mean Diameter - **dp** [:math:`m`]
* Sphericity - **psi** [:math:`-`]
* Minimum Fluidization Voidage - e_mf [:math:`-`]
* Minimum Fluidization Velocity - v_mf [:math:`m/s`]
* Minimum Particle Voidage - **e_min** [:math:`-`]
* Mass Specific Surface Area - **ap** [:math:`m^2/kg`]

Other Variables
^^^^^^^^^^^^^^^

* Elemental composition - **elem_comp** [mol/mol]
* Heat - **Q** [J]
* Work - **W** [J]

Model Initialization
--------------------

All models should include a default initialization routine which is capable of initializing the model from a blank state to a state where it can be successfully solved as part of a flowsheet. Initialization routines should ideally allow the user to specify initial guesses for the state of their system to reduce the amount of effort required to initialize and solve their flowsheets.
