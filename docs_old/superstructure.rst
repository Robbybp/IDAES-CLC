==================================================
Conceptual Design using the IDAES Pyosyn Framework
==================================================

This guide to the IDAES conceptual design framework will illustrate how to develop and optimize models to solve chemical engineering design problems.
This draft of the guide was written based upon commit *d4f6059*, though much of it will also be true for the stable commit *48cc07e*.

To checkout a specific commit, you can simply use ``git checkout 48cc07e``.
This will put you in a detached head mode.
To make a new local branch at that commit, use ``git checkout -b new_branch_name``.

.. contents:: Contents
    :depth: 2

Base Classes
============

A handful of base classes provide necessary and useful functionality when building your models.
The most important of these base classes are listed below:

FlowsheetModel
--------------

Each superstructure synthesis problem should declare a class inheriting from FlowsheetModel that provides information on how to calculate the objective function.

.. module:: idaes_models.core.flowsheet_model 

.. autoclass:: FlowsheetModel

UnitModel
---------

Process units such as reactors, flash drums, compressors, etc. should inherit from UnitModel

.. module:: idaes_models.core.unit_model

# TODO Highlight some key attributes here
equip_exists
equip block
lin_cuts
oa_cuts
units

.. autoclass:: UnitModel

LOA
---

The LOA module provides code for automatically executing the logic-based outer approximation optimization algorithm.

.. automodule:: idaes_models.core.loa
    :members:

Examples
========

Water Treatment Model
---------------------

The water treatment model uses three unit models: Feed, Reactor, and Sink. 
In this example, the objective is to minimize the cost of treating fixed contaminant loads from incoming water streams.
The decision is between one or multiple treatment units in parallel or sequential order.
Here, for convenience, the treatment units are modeled by a unit named *Reactor*.
However, in practice, various forms of filtration, digestion, etc. may be used.
The modeler is free to name their units as they please.

The steps for setting up and solving the water treatment model are thus:

1. Define problem
^^^^^^^^^^^^^^^^^

First, we create the class WaterModel, inheriting from FlowsheetModel.
The objective functions for the problem are defined here.
All modeling objects for the problem are also stored in this class.

.. module:: idaes_models.process.conceptd.water.flowsheet

.. autoclass:: WaterModel
    :members:

2. Create unit models
^^^^^^^^^^^^^^^^^^^^^

Next, we create and define the Feed, Reactor, and Sink unit models.

.. module:: idaes_models.unit.water_net.feed

.. autoclass:: Feed

.. module:: idaes_models.unit.water_net.reactor

.. autoclass:: Reactor

.. module:: idaes_models.unit.water_net.sink

.. autoclass:: Sink

3. Build and connect the model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next step is to build and connect the model. Most of the code for this can be found in the build_model function:

.. module:: idaes_models.process.conceptd.water.flowsheet

.. autofunction:: build_model

The code within this function consists of a section importing data from external sources followed by the creation of the WaterModel, addition of units, and construnction of the units. 
After data import, the first commands are::

  m = WaterModel()
  m.comps = comps
  m.max_flow = 300

This creates a new WaterModel instance, specifies the chemical components that are relevant to the model, and introduces a maximum flowrate for all streams in the flowsheet.

Next, the command ``m.add_unit()`` is used to add the appropriate unit models to the superstructure.
Notice that each unit takes arguments *name* and *parent*.
The name should be unique for all units belonging to the same parent component.
The parent should be the object to which the unit is being added.
It may seem redundant to have ``m.add_unit(Unit(parent=m))``, but it is currently necessary.

.. module:: idaes_models.core.process_base

.. class:: ProcessBase

  .. automethod:: add_unit

Next, connections between units are defined using ``m.connect(from_unit, to_unit)``.
By default, *connect* looks for a port named *outlet* on the from_unit, and a port named *inlet* on the to_unit.
If a unit has multiple possible inlets or outlets, the correct port can be specified using the optional *to_port* or *from_port* arguments:

.. module:: idaes_models.core.process_base

.. class:: ProcessBase

  .. automethod:: connect

The next four commands invoke methods on the unit models to create the optimization model objects::

  m.build_units()
  m.build_links()
  m.expand_connectors()
  propagate_var_fix(m)

The first command executes the ``build()`` command on each unit model.
The next command creates Pyomo Connector objects to link together the unit models.
These Connector objects are then expanded to normal equality constraints by the next command.
Finally, the model is examined for ``a = b`` constraints in which one of the two variables is fixed.
In that case, the other variable linked by equality is also fixed to the same value.

4. Solve the model
^^^^^^^^^^^^^^^^^^

Finally, we solve the conceptual design problem.
The high-level code for doing this can be found in the ``main()`` function:

.. module:: idaes_models.process.conceptd.water.main

.. autofunction:: main

The first command ``m = build_model()`` calls the previously described code.
Next, the code iterates through the superstructure units and applies a linearization strategy::

  for o in itervalues(m.units):
      o.apply_linear_relaxations()
      o.apply_OA_strategy(oa_ports=True)

The first command tells each unit to generate rigorous linear relaxations of nonlinear functions.
The second command tells each unit that OA cuts should be generated for nonlinear equations, including for the interconnection units (ports).

Next, are invocations to the logic-based algorithms::

  do_LOA(m, tol=100)
  do_GLOA(m, tol=100, iterlim=32, contract_bounds=True, do_self_proj=False)

These functions call the respective algorithms in the loa module described above.
Note that in this configuration, the solution to the LOA algorithm is used to initialize the GLOA algorithm, providing an initial upper bound (assuming minimization).
However, either one of the two algorithms can be started and will execute independently of the other.
