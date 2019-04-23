"""
Provides a utility function for comparing two models
"""
from __future__ import division
from pyomo.environ import Block, Var
from math import log10, floor
from six import iteritems
from collections import OrderedDict

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


def compare(m1, m2, **kwargs):
    """
    The idea here is to go block by block through these two models,
    printing out the differences in the variables.
    Some recursion to be expected.
    If m1name and/or m2name is assigned, the corresponding model will be
    renamed to the given string
    """
    m1.name = kwargs.pop('m1name', m1.name)
    m2.name = kwargs.pop('m2name', m2.name)
    diffs = compare_block(m1, m2, **kwargs)
    diff_tups = sorted(iteritems(diffs), key=lambda x: abs(x[1]['diff']))
    if kwargs.pop('display', False) is True:
        for var_id, d_obj in diff_tups:
            print(
                '{}={{{}: {}, {}: {}, diff: {}, rdiff: {}}}'.
                format(var_id,
                       d_obj['name'][0], d_obj['vals'][0],
                       d_obj['name'][1], d_obj['vals'][1],
                       d_obj['diff'], d_obj['rdiff']))
    return OrderedDict(diff_tups)


def compare_block(b1, b2, **kwargs):
    p1names = kwargs.pop('p1names', [])
    p2names = kwargs.pop('p2names', [])
    p1names.append(b1.local_name)
    p2names.append(b2.local_name)

    do_print = kwargs.get('print', False)
    descend_into = kwargs.pop('descend_into', True)
    kwargs.setdefault('diff_objs', {})

    # Get the names of all variables
    b1_varnames = set(b1.component_map(ctype=Var, active=True).iterkeys())
    b2_varnames = set(b2.component_map(ctype=Var, active=True).iterkeys())
    common_vars = b1_varnames & b2_varnames
    if do_print:
        b1_extra_vars = b1_varnames - b2_varnames
        b2_extra_vars = b2_varnames - b1_varnames
        for name in b1_extra_vars:
            print('.'.join(p1names) + ' has extra var ' + name)
        for name in b2_extra_vars:
            print('.'.join(p2names) + ' has extra var ' + name)
    for varname in common_vars:
        compare_var(
            b1.find_component(varname), b2.find_component(varname),
            list(p1names), list(p2names), **kwargs)

    # build a list of subblocks of each model
    b1_subblock_names =\
        set(b1.component_map(ctype=Block, active=True).iterkeys())
    b2_subblock_names =\
        set(b2.component_map(ctype=Block, active=True).iterkeys())
    common_blocks = b1_subblock_names & b2_subblock_names
    if do_print:
        b1_extra_blocks = b1_subblock_names - b2_subblock_names
        b2_extra_blocks = b2_subblock_names - b1_subblock_names
        for name in b1_extra_blocks:
            print('.'.join(p1names) + ' has extra block ' + name)
        for name in b2_extra_blocks:
            print('.'.join(p2names) + ' has extra block ' + name)
    if descend_into:
        for blockname in common_blocks:
            compare_block(
                b1.find_component(blockname), b2.find_component(blockname),
                p1names=list(p1names), p2names=list(p2names), **kwargs)
    return kwargs['diff_objs']

def compare_var(v1, v2, p1names, p2names, **kwargs):
    tol = kwargs.pop('tol', 1E-4)  # tolerance
    sigfigs = kwargs.pop('sigfigs', 3)  # significant figures to display
    # precision to display on relative difference
    rdiff_precision = kwargs.pop('rdiff_precision', 2)
    do_print = kwargs.pop('print', False)
    diff_objs = kwargs.get('diff_objs', {})

    # Define top-level names
    t1name = p1names.pop(0)
    t2name = p2names.pop(0)
    # Variable qualified names (prepending containing blocks)
    v1_qualified_name = '.'.join(p1names + [v1.local_name])
    v2_qualified_name = '.'.join(p2names + [v2.local_name])
    assert(v1_qualified_name == v2_qualified_name)

    common_indx = v1._index & v2._index

    for i in common_indx:
        if v2[i].value is None or v1[i].value is None:
            continue
        diff = v2[i].value - v1[i].value
        rdiff = diff / v1[i].value if v1[i].value != 0 else 0
        if abs(diff) > tol:
            # TODO might want to introduce a relative tolerance too

            v1_rounded = _apply_sigfigs(v1[i].value, sigfigs)
            v2_rounded = _apply_sigfigs(v2[i].value, sigfigs)
            diff_rounded = _apply_sigfigs(diff, sigfigs)
            rdiff_rounded = round(rdiff, rdiff_precision)

            # index string
            istr = '[' + i + ']' if i is not None else ''
            var_id = v1_qualified_name + istr

            # an object describing the difference between the vars
            diff_objs[var_id] = {
                'name': (t1name, t2name),
                'vals': (v1_rounded, v2_rounded),
                'diff': diff_rounded, 'rdiff': rdiff_rounded
            }

            if do_print:
                print(
                    '{}={{{}: {}, {}: {}, diff: {}, rdiff: {}}}'.
                    format(var_id, t1name, v1_rounded, t2name, v2_rounded, diff_rounded, rdiff_rounded))
    # TODO report indices in one variable but not the other
    return diff_objs

def _apply_sigfigs(num, sigfigs):
    if num == 0:
        return 0
    return round(num, -int(floor(log10(abs(num)))) + sigfigs - 1)
