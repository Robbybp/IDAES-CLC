"""
This module provides utility functions for cut generation.
"""
from __future__ import division
from .var import none_if_empty as ne
from pyomo.environ import Var, Param, Block, Constraint, Set
from pyomo.core.base.constraint import SimpleConstraint

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def self_proj_var_rule(main_var, disagg_var, var_set, lbda):
    def var_rule(_, *sets):
        if main_var[ne(sets)].fixed:
            return Constraint.NoConstraint
        else:
            return main_var[ne(sets)] == sum(
                disagg_var[sets, vs] * lbda[vs] for vs in var_set)
    return var_rule

def count_vars(old_block):
    """Count the number of non-fixed variables to clone"""
    count = sum(
        1 for o in old_block.component_data_objects(ctype=Var, active=True, descend_into=False)
        if not o.fixed)
    for old_sub_block in old_block.component_objects(ctype=Block, active=True, descend_into=False):
        if 'equip' not in old_sub_block.name:
            # do not include linear constraints block or ports
            continue
        # clone block
        count += count_vars(old_sub_block)
    return count

def clone_block(old_block, new_block, var_set, lbda):
    """This function acts similarly to the built-in Pyomo clone function,
    but excludes entries that are undesired for this platform."""
    clone_block_vars(old_block, new_block, var_set, lbda)
    clone_block_params(old_block, new_block)
    clone_block_sets(old_block, new_block)
    clone_block_constraints(old_block, new_block, var_set)
    for old_sub_block in old_block.\
            component_objects(ctype=Block, active=True, descend_into=False):
        if 'equip' not in old_sub_block.name:
            # do not include linear constraints block or ports
            continue
        # clone block
        new_sub_block = Block()
        setattr(new_block, old_sub_block.local_name, new_sub_block)
        clone_block(old_sub_block, new_sub_block, var_set, lbda)

def clone_block_vars(old_block, new_block, var_set, lbda):
    for old_var in old_block.\
            component_objects(ctype=Var, active=True, descend_into=False):
        if all(old_var[i].fixed for i in old_var._index):
            # variable is completely fixed. Do not disaggregate; simply copy.
            if old_var.is_indexed():
                new_var = Var(old_var._index)
            else:
                new_var = Var()
            setattr(new_block, old_var.local_name, new_var)
            for i in old_var._index:
                copy_var_data(new_var[i], old_var[i])
        else:
            # create new disaggregated variables
            if old_var.is_indexed():
                new_var = Var(old_var._index, var_set)
                setattr(new_block, old_var.local_name, new_var)
                for i, vs in old_var._index * var_set:
                    copy_var_data(new_var[i, vs], old_var[i])
                # create new aggregated variables
                new_agg_var = Var(old_var._index)
                setattr(new_block, '_' + old_var.local_name + '_clone', new_agg_var)
                for i in old_var._index:
                    copy_var_data(new_agg_var[i], old_var[i])
                # create new aggregate constraints
                setattr(
                    new_block, '_' + old_var.local_name + '_clone_defn',
                    Constraint(old_var._index, rule=self_proj_var_rule(
                        new_agg_var, new_var, var_set, lbda)))
            else:
                # if old_var is a SimpleVar, do not use index i
                new_var = Var(var_set)
                setattr(new_block, old_var.local_name, new_var)
                for vs in var_set:
                    copy_var_data(new_var[vs], old_var)
                # create new aggregated variable
                new_agg_var = Var()
                setattr(new_block, '_' + old_var.local_name + '_clone', new_agg_var)
                copy_var_data(new_agg_var, old_var)
                # create new aggregate constraints
                setattr(
                    new_block, '_' + old_var.local_name + '_clone_defn',
                    Constraint(rule=self_proj_var_rule(
                        new_agg_var, new_var, var_set, lbda)))

def copy_var_data(new_var_data, old_var_data):
    new_var_data.domain = old_var_data.domain
    new_var_data.setlb(old_var_data.lb)
    new_var_data.setub(old_var_data.ub)
    new_var_data.value = old_var_data.value
    new_var_data.fixed = old_var_data.fixed


def clone_block_params(old_block, new_block):
    for old_param in old_block.component_objects(ctype=Param, active=True, descend_into=False):
        # clone param
        new_param = Param(
            old_param._index, initialize=old_param.extract_values())
        setattr(new_block, old_param.local_name, new_param)

def clone_block_sets(old_block, new_block):
    for old_set in old_block.component_objects(ctype=Set, active=True, descend_into=False):
        if old_set.local_name.endswith('_index'):
            continue
        # clone param
        new_param = Set(
            old_set._index, initialize=old_set.data())
        setattr(new_block, old_set.local_name, new_param)

def clone_block_constraints(old_block, new_block, var_set):
    for old_constr in old_block.\
            component_objects(ctype=Constraint, active=True, descend_into=False):
        # clone constraint
        if isinstance(old_constr, SimpleConstraint):
            new_constr = Constraint(var_set, rule=old_constr.rule)
        else:  # old constraint is indexed
            new_constr = Constraint(
                list(old_constr._index), var_set, rule=old_constr.rule)
        setattr(new_block, old_constr.local_name, new_constr)

def get_sum_sq_diff(old_block, new_block):
    sum_sq_diff = 0
    for old_var in old_block.\
            component_objects(ctype=Var, active=True, descend_into=False):
        new_var = getattr(new_block, '_' + old_var.local_name + '_clone', None)
        # new_var = getattr(new_block, old_var.local_name + '_clone', None)
        if new_var is not None:
            # print(old_var.local_name)
            sum_sq_diff += sum(
                (new_var[i] - old_var[i].value) ** 2
                for i in new_var._index
                if new_var is not None and not new_var[i].fixed)
    # recurse through valid subblocks
    for new_sub_block in new_block.\
            component_map(ctype=Block, active=True).itervalues():
        old_sub_block = getattr(old_block, new_sub_block.local_name, None)
        if old_sub_block is not None:
            sum_sq_diff += get_sum_sq_diff(old_sub_block, new_sub_block)
    return sum_sq_diff
