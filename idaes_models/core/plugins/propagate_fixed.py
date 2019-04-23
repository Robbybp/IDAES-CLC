import pyomo.environ as pe

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

# note: see bigm transformation plugin for example of how to make this a plugin

def propagate_var_fix(block, tmp=False):
    """Propagates variable fixing for equalities of type x = y. If x is fixed
    and y is not fixed, then this function will fix y to the value of x.

    If temporary, stores the variable UIDs in a set '_tmp_propagate_fixed'
    attached to the block.

    Args:
        block (TYPE): The block for which to find variables to propagate fixing
        tmp (bool, optional): Whether the variable fixing will be temporary

    Returns:
        None

    Raises:
        RuntimeError: if two fixed variables x = y have different values.
    """
    from pyomo.core.base.expr import _SumExpression
    from pyomo.core.base.var import _GeneralVarData
    if not hasattr(block, '_tmp_propagate_fixed'):
        block._tmp_propagate_fixed = set()
    var_map = {}
    fixed_vars = set()
    eq_var_map = {}
    for constr in block.component_data_objects(ctype=pe.Constraint, active=True, descend_into=True):
        if constr.lower == 0 and constr.upper == 0 \
                and isinstance(constr.body, _SumExpression) \
                and len(constr.body._args) == 2 \
                and 1 in constr.body._coef \
                and -1 in constr.body._coef \
                and isinstance(constr.body._args[0], _GeneralVarData) \
                and isinstance(constr.body._args[1], _GeneralVarData) \
                and constr.body._const == 0:
            v1 = constr.body._args[0]
            v2 = constr.body._args[1]
            if v1.fixed:
                fixed_vars.add(id(v1))
            if v2.fixed:
                fixed_vars.add(id(v2))
            var_map.update({id(v1): v1, id(v2): v2})
            set1 = eq_var_map.get(id(v1), set([id(v1)]))
            set2 = eq_var_map.get(id(v2), set([id(v2)]))
            union = set1.union(set2)
            for vID in union:
                eq_var_map[vID] = union
    processed = set()
    for varID in fixed_vars:
        if varID not in processed:
            var_val = pe.value(var_map[varID])
            eq_set = eq_var_map[varID]
            for v in eq_set:
                if var_map[v].fixed and pe.value(var_map[v]) != var_val:
                    raise RuntimeError('Variables {} and {} have conflicting fixed values of {} and {}, but are linked by equality constraints.'.format(var_map[v1].name, var_map[v2].name, pe.value(var_map[v1]), pe.value(var_map[v2])))
                elif not var_map[v].fixed:
                    var_map[v].fix(var_val)
                    if tmp:
                        block._tmp_propagate_fixed.add(
                            pe.ComponentUID(var_map[v]))
            processed |= eq_set

def reset_propagated_var_fix(block):
    for varUID in block._tmp_propagate_fixed:
        var = varUID.find_component_on(block.model())
        var.unfix()
    block._tmp_propagate_fixed.clear()
