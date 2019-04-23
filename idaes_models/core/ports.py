# -*- coding: utf-8 -*-
from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, Constraint, Block, NonNegativeReals, Set, value, PositiveReals, Param, Binary, Connector, ComponentUID
from idaes_models.core.util.mccormick import add_mccormick_relaxation
from idaes_models.core.util.var import none_if_empty as ne, lb, ub, wrap_var, unwrap_var
from bunch import Bunch
from six import iteritems

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

class UnitPort(UnitModel):
    """Port for a process unit. Acts as either a multi-purpose mixer or
    splitter.

    Attributes:
        F_unit (Var): Pyomo Var describing the total flow from this unit port
            to/from its parent unit
        fc_unit (Var): Pyomo Var describing the component flow from this unit
            port to/from its parent unit
        P_unit (Var): Pyomo Var describing the pressure associated with its
            parent unit
        single_choice (bool): Binary describing whether only one inlet stream may be activated at a time.
        T_unit (Var): Pyomo Var describing the temperature associated with its
            parent unit
        unit_port (Connector): Pyomo Connector linking this unit port with its
            parent unit variables
        vf_unit (Var): Pyomo Var representing the parent unit vapor fraction
            to/from this unit port
        y_unit (Var): Pyomo Var describing the composition of flows between
            this unit port and its parent
    """
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('unit_port')
        self.single_choice = kwargs.pop('single_choice', False)
        unit_vars = ('F', 'T', 'P', 'vf', 'y', 'fc')
        self._unit_vars = Bunch({v: kwargs.pop(v, None) for v in unit_vars})
        kwargs['delay_construct'] = True
        super(UnitPort, self).__init__(*args, **kwargs)

    def __init2__(self, *args, **kwargs):
        super(UnitPort, self).__init2__(*args, **kwargs)
        self.build_unit_facing_vars()

    def build_unit_facing_vars(self):
        if self._unit_vars.P is not None:
            self.P_unit = Var(domain=NonNegativeReals, initialize=1.0, bounds=(0, self.max_P))
        if self._unit_vars.T is not None:
            self.T_unit = Var(domain=PositiveReals, initialize=300, bounds=(self.min_T, self.max_T))
        self.fc_unit = Var(self.comps, domain=NonNegativeReals, initialize=1.0, bounds=(0, self._max_flow))
        self.F_unit = Var(domain=NonNegativeReals, initialize=1.0, bounds=(0, self._max_flow))
        self.y_unit = Var(self.comps, domain=NonNegativeReals, initialize=0.5, bounds=(0, 1), doc="Mole fraction to/from the unit")
        if self._unit_vars.vf is not None:
            self.vf_unit = Var(domain=NonNegativeReals, initialize=0.5, bounds=(0, 1))
        if all(v is not None for v in (
                self._unit_vars.fc, self._unit_vars.F, self._unit_vars.y)):
            # indeterminant flow scheme. Use default
            self._flow_scheme = 'fc'
        elif self._unit_vars.fc is not None:
            self._flow_scheme = 'fc'
        elif self._unit_vars.F is not None and self._unit_vars.y is not None:
            self._flow_scheme = 'Fy'
        else:
            raise Exception('Unknown flow scheme. Need to specify F and y, or fc')
        unit_c = Connector()
        port_c = Connector()
        for k, v in iteritems(self._unit_vars):
            if v is not None:
                unit_c.add(v, name=k)
                port_c.add(getattr(self, k + '_unit'), name=k)
        # attach the new connector to the parent block
        setattr(self.parent_unit._obj, '{}_port'.format(self.unit_name), unit_c)
        # attach the new connector to the UnitPort block
        self.unit_port = port_c
        # build the link between the parent and this UnitPort
        link_name = "var_link_with_{}".format(self.unit_name)
        self.parent_unit.link_map[link_name] = {'var_port': unit_c, 'unit_port': port_c}

    def generate_cut_gen_problem(self):
        # Do nothing for now
        pass

    def fathom(self):
        """Overrides superclass fathom to make sure that parent UnitModel also registers the fathom. Needs to have logic ensuring that an infinite recursion loop is not formed.

        Returns:
            None
        """
        if not self._fathomed:
            self._fathomed = True
            self.parent_unit.fathom()
        else:
            super(UnitPort, self).fathom()

    def reset_introspect_fixed(self):
        """Resets variables fixed and constraints deactivated during
        introspection

        Returns:
            None
        """
        for varUID in self._tmp_introspect_fixed:
            var = varUID.find_component_on(self.model())
            var.unfix()
        self._tmp_introspect_fixed.clear()
        for conUID in self._tmp_introspect_deactivated:
            condata = conUID.find_component_on(self.model())
            condata.activate()


@ProcBlock("InPort")
class _InPort(UnitPort):
    """
    Inlet port to a process unit. Acts as a multi-purpose mixer.

    Attributes:
        Cp (dict): Dictionary of heat capacities for each component. By default, all components are assigned a Cp of 1.
        FT (Var): Pyomo Var describing the product of total flow and temperature
        in_streams (Set): Pyomo Set of inlet streams
        in_streams_list (list): List of strings corresponding to inlet streams.
            This list is usually populated by calling the add_source function.
        lin_cuts (Block): see UnitModel attribute description
        lin_cuts_nsegs (int): number of piecewise segments to use for
            constructing linear relaxations
        out_stream (Set): Pyomo Set containing the name of the outlet stream
            from this inlet port.
        stream_active (Var): Boolean variable describing if a given inlet
            stream is active
    """
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('in_port')
        self.in_streams_list = kwargs.pop('in_stream_names', [])
        self._Cp = kwargs.pop('Cp', {})
        super(_InPort, self).__init__(*args, **kwargs)

    def add_source(self, from_port):
        self.in_streams_list.append(from_port.parent_unit.unit_name)

    def build(self):
        if len(self.in_streams_list) >= 1:
            self.base_stream = self.in_streams_list[0]

        self.equip_exists_var = wrap_var(getattr(self.parent_unit, 'equip_exists', 1.0))

        super(_InPort, self).build()

    @property
    def equip_exists(self):
        return value(self.equip_exists_var)

    def _build_unit_sets(self):
        super(_InPort, self)._build_unit_sets()
        self.in_streams = Set(initialize=self.in_streams_list)
        self.out_stream = Set(initialize=['unit'])

    def _build_unit_params(self):
        super(_InPort, self)._build_unit_params()
        if not self._Cp:
            self._Cp = {s: 1.0 for s in self.in_streams}
            self._Cp['unit'] = 1.0
        self.Cp = Param(self.in_streams | self.out_stream, initialize=self._Cp)

    def _build_unit_vars(self):
        super(_InPort, self)._build_unit_vars()

        for stream in sorted(self.in_streams):
            if self._unit_vars.P is not None:
                P_var = Var(domain=NonNegativeReals, initialize=1.0, bounds=(0, self.max_P))
                setattr(self, 'P_' + stream, P_var)

            if self._unit_vars.T is not None:
                T_var = Var(domain=PositiveReals, initialize=300, bounds=(self.min_T, self.max_T))
                setattr(self, 'T_' + stream, T_var)

            fc_var = Var(self.comps, domain=NonNegativeReals, initialize=1.0, bounds=(0, self._max_flow))
            setattr(self, 'fc_' + stream, fc_var)

            F_var = Var(domain=NonNegativeReals, initialize=1.0, bounds=(0, self._max_flow))
            setattr(self, 'F_' + stream, F_var)

            y_var = Var(self.comps, domain=NonNegativeReals, initialize=0.5, bounds=(0, 1))
            setattr(self, 'y_' + stream, y_var)

            vf_var = Var(domain=NonNegativeReals, initialize=0.5, bounds=(0, 1))
            setattr(self, 'vf_' + stream, vf_var)

        if self._unit_vars.T is not None and len(self.in_streams) >= 2:
            self.FT = Var(self.in_streams | self.out_stream, domain=NonNegativeReals, initialize=0)

        if self._unit_vars.T is not None and len(self.in_streams) >= 2 and self.single_choice:
            self.stream_active = Var(self.in_streams, domain=Binary)

        # Initialize variables
        # TODO may want to reassess how to initialize
        for stream in self.in_streams:
            getattr(self, 'F_' + stream).set_value(sum(value(getattr(self, 'fc_' + stream)[c]) for c in self.comps))

    def _build_unit_constraints(self):
        super(_InPort, self)._build_unit_constraints()
        b = self
        equip = self.equip

        @equip.Constraint(self.comps)
        def fc_mass_balance(equip, c):
            """f_unit,c = Σ_s( f_s,c )"""
            return b.fc_unit[c] == sum(getattr(b, 'fc_' + s)[c] for s in self.in_streams)
        if len(self.in_streams) < 1:
            # Only apply mass balance if there are inlet streams
            equip.fc_mass_balance.deactivate()

        @equip.Constraint()
        def F_mass_balance(equip):
            """F_unit = Σ_s( F_s )"""
            return b.F_unit == sum(getattr(b, 'F_' + s) for s in self.in_streams)
        if not self._flow_scheme == 'Fy' or len(self.in_streams) < 1:
            equip.F_mass_balance.deactivate()

        @equip.Constraint(self.in_streams | self.out_stream)
        def total_flow(equip, s):
            """F_unit = Σ_c( f_s,c )"""
            return getattr(b, 'F_' + s) == sum(getattr(b, 'fc_' + s)[c] for c in self.comps)
        if not self._flow_scheme == 'fc':
            equip.total_flow.deactivate()

        if self._flow_scheme == 'Fy':
            @equip.Constraint(self.in_streams)
            def sum_y(equip, s):
                """Σ_c( y_s,c ) = 1"""
                return sum(getattr(b, 'y_' + s)[c] for c in self.comps) == 1

        @equip.Constraint(self.in_streams | self.out_stream, self.comps)
        def flow_comps(equip, s, c):
            """fc_s,c = F_s × y_s,c
            or fc_unit_c = F_unit × y_unit,c"""
            return getattr(b, 'fc_' + s)[c] == getattr(b, 'F_' + s) * getattr(b, 'y_' + s)[c]
        if len(self.in_streams) < 2:  # proper
            for s, c in self.in_streams * self.comps:
                equip.flow_comps[s, c].deactivate()

        if len(self.in_streams) == 1:
            @equip.Constraint(self.comps)
            def y_eq(equip, c):
                return getattr(b, 'y_' + b.base_stream)[c] == b.y_unit[c]

        if self._flow_scheme == 'fc':
            # If the flow scheme is component flows, do not calculate Fy
            # quantities unless necessary.
            equip.flow_comps.deactivate()
            try:
                equip.y_eq.deactivate()
            except AttributeError:
                pass

        if self._unit_vars.P is not None:
            @equip.Constraint(self.in_streams)
            def P_equal(equip, s):
                return getattr(b, 'P_' + s) == b.P_unit

        if self._unit_vars.T is not None and len(self.in_streams) >= 2:
            @equip.Constraint()
            def energy_balance(equip):
                return b.F_unit * b.Cp['unit'] * b.T_unit == sum(getattr(self, 'F_' + s) * b.Cp[s] * getattr(self, 'T_' + s) for s in b.in_streams)

        if self._unit_vars.T is not None and len(self.in_streams) == 1:
            @equip.Constraint()
            def T_equal(equip):
                return b.T_unit == getattr(self, 'T_' + self.base_stream)

        # Variable connectors
        for stream in self.in_streams:
            port_c = Connector()
            for k, v in iteritems(self._unit_vars):
                if v is not None:
                    port_c.add(getattr(self, k + '_' + stream), name=k)
            setattr(self, stream + '_port', port_c)

    def get_external_link(self, link):
        return getattr(self, link['from'].unit_name + '_port')

    def get_flow_vars(self):
        yield self.fc_unit
        yield self.F_unit
        for s in self.in_streams:
            yield getattr(self, 'fc_' + s)
            yield getattr(self, 'F_' + s)

    def introspect_flows(self):
        """This function examines the flow values of attached connections for
        those that are deactivated, and fixes the corresponding flow-related
        variables to zero. These variables are stored via UID in a temporary
        set '_tmp_introspect_fixed' so that they can be unfixed later.

        This function also deactivates certain constraints that are no longer
        needed in the case of a deactivated flow. These constraints are stored
        via UID in a temporary set '_tmp_introspect_deactivated' so that they
        can be reactivated later.

        Returns:
            None
        """
        if not hasattr(self, '_tmp_introspect_fixed'):
            self._tmp_introspect_fixed = set()
        if not hasattr(self, '_tmp_introspect_deactivated'):
            self._tmp_introspect_deactivated = set()
        for s in self.in_streams:
            # if all the component flows are fixed to zero, then the stream is deactivated.
            fc = getattr(self, 'fc_' + s)
            F = getattr(self, 'F_' + s)
            if all(value(fc[c]) == 0 and fc[c].fixed for c in self.comps) or \
                    (value(F) == 0 and F.fixed):
                # Deactivate associated variables
                for c in self.comps:
                    if not fc[c].fixed:
                        self._tmp_introspect_fixed.add(ComponentUID(fc[c]))
                        fc[c].fix(0)
                if not F.fixed:
                    self._tmp_introspect_fixed.add(ComponentUID(F))
                    F.fix(0)
                # If there are equations to deactivate, do it here.
                # for c in self.comps:
                #     try:
                #         condata = self.equip.exit_flows[s, c]
                #     except AttributeError:
                #         break
                #     except KeyError:
                #         # Maybe if there's a case that a certain component
                #         # won't be defined, this code would be a problem?
                #         break
                #     if condata.active:
                #         condata.deactivate()
                #         self._tmp_introspect_deactivated.add(ComponentUID(condata))

    def apply_linear_relaxations(self, nsegs=1):
        from idaes_models.core.util.mccormick import add_mccormick_relaxation
        from idaes_models.core.util.var import tighten_mc_var
        from functools import partial
        tighten_mc = partial(tighten_mc_var, block_bounds=self.block_bounds)

        b = self
        lin_cuts = self.lin_cuts
        self.lin_cuts_nsegs = nsegs

        if len(self.in_streams) >= 1:
            lin_cuts.mc_flow_unit = Block()
            for c in self.comps:
                tighten_mc(b.fc_unit[c], b.F_unit, b.y_unit[c])
                add_mccormick_relaxation(lin_cuts.mc_flow_unit, b.fc_unit[c], b.F_unit, b.y_unit[c], nsegs, c, 1.0, block_bounds=self.block_bounds)

        if len(self.in_streams) >= 2 and not self.single_choice:
            lin_cuts.mc_flow_comps = Block()
            for s, c in self.in_streams * self.comps:
                tighten_mc(getattr(b, 'fc_' + s)[c], getattr(b, 'F_' + s), getattr(b, 'y_' + s)[c])
                add_mccormick_relaxation(lin_cuts.mc_flow_comps, getattr(b, 'fc_' + s)[c], getattr(b, 'F_' + s), getattr(b, 'y_' + s)[c], nsegs, (s, c), 1.0, block_bounds=self.block_bounds)
        elif len(self.in_streams) >= 2 and self.single_choice:
            @lin_cuts.Constraint(self.in_streams, self.comps)
            def fc_lb(lin_cuts, s, c):
                return getattr(b, 'fc_' + s)[c] >= b.fc_unit[c] - (ub(b.fc_unit[c]) - lb(getattr(b, 'fc_' + s)[c])) * (1 - b.stream_active[s])

            @lin_cuts.Constraint(self.in_streams, self.comps)
            def fc_ub(lin_cuts, s, c):
                return getattr(b, 'fc_' + s)[c] <= ub(b.fc_unit[c]) * b.stream_active[s]

            if self._flow_scheme == 'Fy':
                # also bound F and y
                @lin_cuts.Constraint(self.in_streams)
                def F_lb(lin_cuts, s):
                    return getattr(b, 'F_' + s) >= b.F_unit - (ub(b.F_unit) - lb(getattr(b, 'F_' + s)) * (1 - b.stream_active[s]))

                @lin_cuts.Constraint(self.in_streams)
                def F_ub(lin_cuts, s):
                    return getattr(b, 'F_' + s) <= ub(b.F_unit) * b.stream_active[s]

                @lin_cuts.Constraint(self.in_streams, self.comps)
                def y_lb(lin_cuts, s, c):
                    return getattr(b, 'y_' + s)[c] >= b.y_unit - (ub(b.y_unit) - lb(getattr(b, 'y_' + s))) * (1 - b.stream_active[s])

                @lin_cuts.Constraint(self.in_streams, self.comps)
                def y_ub(lin_cuts, s, c):
                    return getattr(b, 'y_' + s)[c] <= ub(b.y_unit) * b.stream_active[s]

        if self._unit_vars.T is not None and len(self.in_streams) >= 2 and not self.single_choice:
            lin_cuts.mc_energy_balance = Block()
            for s in self.in_streams | self.out_stream:
                tighten_mc(self.FT[s], getattr(self, 'F_' + s), getattr(self, 'T_' + s))
                add_mccormick_relaxation(lin_cuts.mc_energy_balance, self.FT[s], getattr(self, 'F_' + s), getattr(self, 'T_' + s), nsegs, s, 1.0, block_bounds=self.block_bounds)

            lin_cuts.linear_energy_balance = Constraint(expr=self.FT['unit'] * self.Cp['unit'] == sum(self.FT[s] * self.Cp[s] for s in self.in_streams))
        elif self._unit_vars.T is not None and len(self.in_streams) >= 2 and self.single_choice:
            def Tout_lb(lin_cuts, s):
                b = lin_cuts.parent_block()
                Tstr = getattr(self, 'T_' + s)
                return b.T_unit >= Tstr - (ub(Tstr) - lb(b.T_unit)) * (1 - b.stream_active[s])
            lin_cuts.Tout_lb = Constraint(self.in_streams, rule=Tout_lb)

            def Tout_ub(lin_cuts, s):
                b = lin_cuts.parent_block()
                Tstr = getattr(self, 'T_' + s)
                return b.T_unit <= Tstr + (ub(b.T_unit) - lb(Tstr)) * (1 - b.stream_active[s])
            lin_cuts.Tout_ub = Constraint(self.in_streams, rule=Tout_ub)

            def flow_ub(lin_cuts, s, c):
                b = lin_cuts.parent_block()
                return getattr(self, 'fc_' + s)[c] <= ub(getattr(self, 'fc_' + s)[c]) * b.stream_active[s]
            lin_cuts.flow_ub = Constraint(self.in_streams, self.comps, rule=flow_ub)

            lin_cuts.single_choice = Constraint(expr=sum(self.stream_active[s] for s in self.in_streams) == 1)

    def reconstruct_envelopes(self):
        self.del_component('lin_cuts')
        self.lin_cuts = Block()
        self.apply_linear_relaxations(nsegs=self.lin_cuts_nsegs)

@ProcBlock("OutPort")
class _OutPort(UnitPort):
    """
    Outlet port from a process unit. Acts as a multi-purpose splitter.

    Attributes:
        base_stream (string): The identifier for the stream excluded from the exit_flows equation to make the system of equations more square. Arbitrarily selected.
        equip_exists_var (Var or float): Pyomo variable corresponding to whether the parent unit exists. If parent unit does not have this variable, then this will be the float 1.0
        in_stream (Set): The inlet stream to this outlet port
        lin_cuts (Block): see description in UnitModel
        lin_cuts_nsegs (int): number of piecewise segments to use for
            constructing linear relaxations
        out_streams (Set): Pyomo Set of outlet streams
        out_streams_list (list): List of string corresponding to outlet
            streams. Add to this list using the add_destination function.
    """
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('out_port')
        self.out_streams_list = kwargs.pop('out_streams', [])
        super(_OutPort, self).__init__(*args, **kwargs)
        # TODO handle in a more elegant way the case of having 1 or 0 streams
        # connected to out port

    def add_destination(self, to_port):
        self.out_streams_list.append(to_port.parent_unit.unit_name)

    def get_external_link(self, link):
        return getattr(self, link['to'].unit_name + '_port')

    def _build_unit_sets(self):
        super(_OutPort, self)._build_unit_sets()
        self.out_streams = Set(initialize=self.out_streams_list)
        self.in_stream = Set(initialize=['unit'])

    def _build_unit_params(self):
        super(_OutPort, self)._build_unit_params()

    def _build_unit_vars(self):
        super(_OutPort, self)._build_unit_vars()
        b = self
        for stream in b.out_streams:
            if self._unit_vars.T is not None:
                T_var = Var(domain=PositiveReals, initialize=300, bounds=(self.min_T, self.max_T))
                setattr(self, 'T_' + stream, T_var)

            if self._unit_vars.P is not None:
                P_var = Var(domain=NonNegativeReals, initialize=1.0, bounds=(0, self.max_P))
                setattr(self, 'P_' + stream, P_var)

            fc_var = Var(self.comps, domain=NonNegativeReals, bounds=(0, self._max_flow), initialize=1.0)
            setattr(self, 'fc_' + stream, fc_var)

            F_var = Var(domain=NonNegativeReals, initialize=1.0, bounds=(0, self._max_flow))
            setattr(self, 'F_' + stream, F_var)

            y_var = Var(self.comps, domain=NonNegativeReals, initialize=0.5, bounds=(0, 1))
            setattr(self, 'y_' + stream, y_var)

            vf_var = Var(domain=NonNegativeReals, initialize=0.5, bounds=(0, 1))
            setattr(self, 'vf_' + stream, vf_var)

        if len(b.out_streams) >= 1:
            b.frac = Var(b.out_streams, domain=NonNegativeReals, initialize=0.5, bounds=(0, 1), doc="Fraction of flow to each outlet stream")

        if self.single_choice and len(b.out_streams) >= 2:
            b.stream_active = Var(b.out_streams, domain=Binary)

        # initialize fc_unit
        for c in self.comps:
            b.fc_unit[c] = sum(value(getattr(self, 'fc_' + s)[c]) for s in b.out_streams)

        b.F_unit = sum(value(b.fc_unit[c]) for c in self.comps)
        # initialize frac
        for s in b.out_streams:
            init_frac = value(getattr(self, 'F_' + s)) / value(b.F_unit) if value(b.F_unit) > 0 else 0.5
            init_frac = min(max(0, init_frac), 1)
            b.frac[s] = init_frac

    def _build_unit_constraints(self):
        super(_OutPort, self)._build_unit_constraints()
        b = self
        equip = b.equip

        # We actually do not enforce this constraint for the nonlinear case,
        # because it causes redundancy in the equations.
        #
        # @equip.Constraint(self.comps)
        # def mass_balance(equip, c):
        #     """f_unit,c = Σ_s( f_s,c )"""
        #     return sum(getattr(b, 'fc_' + s)[c] for s in b.out_streams) == b.fc_unit[c]

        @equip.Constraint()
        def F_mass_balance(equip):
            return sum(getattr(b, 'F_' + s) for s in b.out_streams) == b.F_unit

        if not ((len(b.out_streams) >= 1 and self._flow_scheme == 'Fy') or
                (len(b.out_streams) == 1 and self._flow_scheme == 'fc')):
            equip.F_mass_balance.deactivate()

        if len(b.out_streams) >= 1:
            @equip.Constraint(b.out_streams, b.comps)
            def exit_flows(equip, s, c, *vs):
                return getattr(b, 'fc_' + s)[c, vs] == b.fc_unit[c, vs] * b.frac[s, vs]
            if self._flow_scheme == 'Fy':
                equip.exit_flows.deactivate()

        if self._flow_scheme == 'Fy' or \
                (self._flow_scheme == 'fc' and len(b.out_streams) == 1):
            @equip.Constraint(self.out_streams, self.comps)
            def y_equal(equip, s, c):
                """y_unit,c = y_s,c"""
                return b.y_unit[c] == getattr(b, 'y_' + s)[c]

            # only enforce if stream is active?

        if self._flow_scheme == 'Fy':
            @equip.Constraint(b.in_stream | b.out_streams)
            def y_sum(equip, s):
                return sum(getattr(b, 'y_' + s)[c] for c in self.comps) == 1

        @equip.Constraint(b.in_stream | b.out_streams)
        def total_flow_defn(equip, s):
            return getattr(b, 'F_' + s) == sum(getattr(b, 'fc_' + s)[c] for c in self.comps)
        if not (self._flow_scheme == 'fc' and len(b.out_streams) >= 2):
            for s in b.out_streams:
                equip.total_flow_defn[s].deactivate()
        if not self._flow_scheme == 'fc':
            equip.total_flow_defn['unit'].deactivate()

        if len(b.out_streams) >= 2:
            @equip.Constraint()
            def frac_sum(equip):
                return sum(b.frac[s] for s in b.out_streams) == 1
        elif len(b.out_streams) == 1:
            b.frac[b.base_stream].fix(1)

        @equip.Constraint(self.out_streams | self.in_stream, self.comps)
        def flow_comps(equip, s, c):
            """fc_s,c = F_s × y_s,c
            or fc_unit_c = F_unit × y_unit,c"""
            return getattr(b, 'fc_' + s)[c] == getattr(b, 'F_' + s) * getattr(b, 'y_' + s)[c]
        if not self._flow_scheme == 'fc':
            for c in self.comps:
                equip.flow_comps['unit', c].deactivate()
        if not ((len(self.out_streams) >= 2 and self._flow_scheme == 'fc') or
                (len(self.out_streams) >= 1 and self._flow_scheme == 'Fy')):
            for s, c in self.out_streams * self.comps:
                equip.flow_comps[s, c].deactivate()

        if self._flow_scheme == 'fc':
            # If flow scheme is component flows, do not calculate Fy quantities
            # unless necessary.
            equip.flow_comps.deactivate()

        if len(b.out_streams) >= 1 and self._unit_vars.T is not None:
            @equip.Constraint(b.out_streams)
            def T_equal(equip, s, *vs):
                b = equip.parent_block()
                return getattr(b, 'T_' + s)[ne(vs)] == b.T_unit

        if len(b.out_streams) >= 1 and self._unit_vars.P is not None:
            @equip.Constraint(b.out_streams)
            def P_equal(equip, s, *vs):
                b = equip.parent_block()
                return getattr(b, 'P_' + s)[ne(vs)] == b.P_unit

        if self._unit_vars.vf is not None:
            @equip.Constraint(b.out_streams)
            def vf_equal(equip, s):
                return getattr(b, 'vf_' + s) == b.vf_unit

        for stream in b.out_streams:
            port_c = Connector()
            for k, v in iteritems(self._unit_vars):
                if v is not None:
                    port_c.add(getattr(self, k + '_' + stream), name=k)
            setattr(self, stream + '_port', port_c)

    def build(self):
        if len(self.out_streams_list) >= 1:
            self.base_stream = self.out_streams_list[0]

        self.equip_exists_var = wrap_var(getattr(self.parent_unit, 'equip_exists', 1.0))

        super(_OutPort, self).build()

    @property
    def equip_exists(self):
        return value(self.equip_exists_var)

    def introspect_flows(self):
        """This function examines the flow values of attached connections for
        those that are deactivated, and fixes the corresponding flow-related
        variables to zero. These variables are stored via UID in a temporary
        set '_tmp_introspect_fixed' so that they can be unfixed later.

        This function also deactivates certain constraints that are no longer
        needed in the case of a deactivated flow. These constraints are stored
        via UID in a temporary set '_tmp_introspect_deactivated' so that they
        can be reactivated later.

        Returns:
            None
        """
        if not hasattr(self, '_tmp_introspect_fixed'):
            self._tmp_introspect_fixed = set()
        if not hasattr(self, '_tmp_introspect_deactivated'):
            self._tmp_introspect_deactivated = set()
        for s in self.out_streams:
            # if all the component flows are fixed to zero, then the stream is deactivated.
            fc = getattr(self, 'fc_' + s)
            F = getattr(self, 'F_' + s)
            xi = self.frac[s]
            if all(value(fc[c]) == 0 and fc[c].fixed for c in self.comps) or \
                    (value(F) == 0 and F.fixed) or \
                    (value(xi) == 0 and xi.fixed):
                # Fix associated variables
                for c in self.comps:
                    if not fc[c].fixed:
                        self._tmp_introspect_fixed.add(ComponentUID(fc[c]))
                        fc[c].fix(0)
                if not F.fixed:
                    self._tmp_introspect_fixed.add(ComponentUID(F))
                    F.fix(0)
                if not xi.fixed:
                    self._tmp_introspect_fixed.add(ComponentUID(xi))
                    xi.fix(0)
                # Deactivate associated constraints
                for c in self.comps:
                    try:
                        condata = self.equip.flow_comps_out[s, c]
                    except AttributeError:
                        break
                    if condata.active:
                        condata.deactivate()
                        self._tmp_introspect_deactivated.add(ComponentUID(condata))
                for c in self.comps:
                    try:
                        condata = self.equip.exit_flows[s, c]
                    except AttributeError:
                        break
                    except KeyError:
                        # Maybe if there's a case that a certain component
                        # won't be defined, this code would be a problem?
                        break
                    if condata.active:
                        condata.deactivate()
                        self._tmp_introspect_deactivated.add(ComponentUID(condata))

    def apply_linear_relaxations(self, nsegs=1):
        b = self
        lin_cuts = self.lin_cuts
        self.lin_cuts_nsegs = nsegs

        def mass_balance(lin_cuts, c, *vs):
            """f_unit,c = Σ_s( f_s,c )"""
            return sum(getattr(b, 'fc_' + s)[c, vs] for s in b.out_streams) == b.fc_unit[c, vs]

        if len(self.out_streams) >= 1:
            lin_cuts.mass_balance = Constraint(self.comps, rule=mass_balance)

        if not self.single_choice and len(self.out_streams) >= 2:
            mc_block = lin_cuts.mc_exit_flows = Block()
            for (c, s) in self.comps * self.out_streams:
                add_mccormick_relaxation(mc_block, getattr(self, 'fc_' + s)[c], self.fc_unit[c], self.frac[s], nsegs, (c, s), unwrap_var(self.equip_exists_var), self.block_bounds)

        if self.single_choice and len(self.out_streams) >= 2:
            def glover_flow_ub(lin_cuts, s, c):
                return getattr(self, 'fc_' + s)[c] <= ub(self.fc_unit[c]) * self.stream_active[s]
            lin_cuts.glover_flow_ub = Constraint(self.out_streams, self.comps, rule=glover_flow_ub)

            def glover_flow_equal_1(lin_cuts, s, c):
                return getattr(self, 'fc_' + s)[c] >= self.fc_unit[c] - (1 - self.stream_active[s]) * ub(self.fc_unit[c])
            lin_cuts.glover_flow_equal_1 = Constraint(self.out_streams, self.comps, rule=glover_flow_equal_1)

            def glover_flow_equal_2(lin_cuts, s, c):
                return getattr(self, 'fc_' + s)[c] <= self.fc_unit[c] - (1 - self.stream_active[s]) * lb(self.fc_unit[c])
            lin_cuts.glover_flow_equal_2 = Constraint(self.out_streams, self.comps, rule=glover_flow_equal_2)

    def get_vars_to_bound(self):
        yield self.fc_unit
        yield self.frac

    def reconstruct_envelopes(self):
        self.del_component('lin_cuts')
        self.lin_cuts = Block()
        self.apply_linear_relaxations(nsegs=self.lin_cuts_nsegs)

    def get_flow_vars(self):
        yield self.fc_unit
        yield self.F_unit
        for s in self.out_streams:
            yield getattr(self, 'fc_' + s)
            yield getattr(self, 'F_' + s)
