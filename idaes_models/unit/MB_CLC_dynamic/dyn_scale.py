from pyomo.environ import (Var, SimpleVar, Expression, Constraint, Suffix)
import pdb

def replace_variables(m, m_ss):
    for var in m.component_objects(Var, descend_into=True):
        # maybe make a function that acts on variables
        # find variable's value at initial and final steady state
        # if variable indexed by time and changes with time
        #    create reference parameter/(local) suffix as final steady state value
        #    create scale export suffix as initial-final ss value
        #    create difference Expression
        # else
        #    create scale export suffix as (initial) value of var

        # maybe a dictionary of var data names (or ids) to Expressions
        print(var.name)

    for e in m.component_data_objects(Constraint, descend_into=True):
        print(e.name)

def is_indexed_by(var, s):
    # returns boolean saying if var is indexed by set s
    if isinstance(var, SimpleVar):
        return False
    n = var.index_set().dimen
    if n == 1:
        if var.index_set() == s:
            return True
        return False
    if n >= 2:
        if s in var.index_set().set_tuple:
            return True
        return False

def get_non_time_index(var, index, t):
    if not is_indexed_by(var, t):
        return index
    if var.index_set() == t:
        return None
    n = var.index_set().dimen
    if n < 2:
        raise ValueError('''Inconsistency getting non-time index.
                            Expected variable to be indexed by multiple sets.''')
    if n >= 2:
        time_loc = -1
        for i in range(0,n):
            if t == var.index_set().set_tuple[i]:
                time_loc = i
        if time_loc == -1:
            raise ValueError('''Inconsistency getting non-time index.
                                Expected to find time in the set_tuple''')
        index_list = []
        for i in range(0,n):
            if i != time_loc:
                index_list.append(index[i])
        index_tuple = tuple(index_list)
        if len(index_tuple) == 1:
            return index_tuple[0]
        elif len(index_tuple) > 1:
            return index_tuple

def create_suffixes(fs):
    m = fs.MB_fuel
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)
    m.ss_init = Suffix(direction=Suffix.LOCAL)
    m.ss_fin = Suffix(direction=Suffix.LOCAL)

def create_scale_values(var, fs, fs_ss_init, fs_ss_fin):
    m = fs.MB_fuel
    m_ss_init = fs_ss_init.MB_fuel
    m_ss_fin = fs_ss_fin.MB_fuel
    time = m.t
    varname = var.local_name

    var_ss_init = eval('m_ss_init.'+varname)
    var_ss_fin = eval('m_ss_fin.'+varname)

    if not is_indexed_by(var, time):
        return

    for index in var:
        nt_index = get_non_time_index(var, index, time)
        m.ss_init[var[index]] = var_ss_init[nt_index].value
        m.ss_fin[var[index]] = var_ss_fin[nt_index].value

        if var_ss_init[nt_index].value == var_ss_fin[nt_index].value:
            m.scaling_factor[var[index]] = var_ss_init[nt_index].value

        elif var_ss_init[nt_index].value != var_ss_fin[nt_index].value:
            m.scaling_factor[var[index]] = (var_ss_init[nt_index].value -
                                            var_ss_fin[nt_index].value)



