from pyomo.environ import (Var, SimpleVar, Expression, Constraint, Suffix,
        Param)
from pyomo.core.expr.current import (SimpleExpressionVisitor, 
        ExpressionReplacementVisitor, ExpressionBase)
from pyomo.core.expr.numvalue import (native_numeric_types,
        NumericValue, NumericConstant)
from pyomo.dae import DerivativeVar
import pdb

class DeviationVisitor(ExpressionReplacementVisitor):
    
    def __init__(self):
        super(DeviationVisitor, self).__init__()

    def visiting_potential_leaf(self, node):
        print('node:', node)
        if node.__class__ in native_numeric_types:
            return True, node

        if node.is_variable_type():
            try: 
                # okay for now, but node might not have parent component but still have dev exp
                # ^ this doesn't seem to make sense, but I have no proof
                if node.parent_component().has_dev_exp == True:
                    print('replacing')
                    return True, node.dev_exp[node.index().expr]
            except AttributeError:
                pass
                
        return False, node

class TestWalker(SimpleExpressionVisitor):

    def __init__(self):
        super(SimpleExpressionVisitor).__init__()

    def visit(self, node):
        print('type:', type(node))
        print(isinstance(node, ExpressionBase))
        #if not type(node) in native_numeric_types and not isinstance(node, NumericConstant):
        #    print('nargs:', node.nargs())
        print(node, '\n')

def walk_tree(expr):
    visitor = TestWalker()
    return visitor.xbfs(expr)

def replace_variables(expr):
    visitor = DeviationVisitor()
    return visitor.dfs_postorder_stack(expr)

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

    if isinstance(var, DerivativeVar):
        return

    if var.index_set().dimen == 1:
        index_sets = [var.index_set()]
    elif var.index_set().dimen >= 2:
        index_sets = list(var.index_set().set_tuple)

    # need to remember which var data objects have deviation variables
    #var.dev = Var(*index_sets,initialize=

    var.has_dev_exp = False

    for index in var:
        nt_index = get_non_time_index(var, index, time)
        m.ss_init[var[index]] = var_ss_init[nt_index].value
        m.ss_fin[var[index]] = var_ss_fin[nt_index].value

        if var_ss_init[nt_index].value == var_ss_fin[nt_index].value:
            m.scaling_factor[var[index]] = var_ss_init[nt_index].value

        elif var_ss_init[nt_index].value != var_ss_fin[nt_index].value:
            m.scaling_factor[var[index]] = (var_ss_init[nt_index].value -
                                            var_ss_fin[nt_index].value)
            var.has_dev_exp = True
            
    def dev_init_rule(m, *args):
        # can a variable indexed by a single set accept a tuple as a subscript?
        # I hope so...
        return var[tuple(args)].value - m.ss_fin[var[tuple(args)]]
    var.dev = Var(*index_sets, initialize=dev_init_rule)

    m.add_component(varname + '_dev', var.dev)
    var_dev = eval('m.'+varname+'_dev')

    for index in var_dev:
        if var_dev[index].value != 0:
            m.scaling_factor[var_dev[index]] = m.scaling_factor[var[index]]

    def dev_exp_rule(m, *args):
        if m.ss_init[var[tuple(args)]] == m.ss_fin[var[tuple(args)]]:
            # not super clear what I should do in this case...
            # returns the original variable, will replace it by itself...
            return var[tuple(args)]
        else:
            return var_dev[tuple(args)] + m.ss_fin[var[tuple(args)]]
                #(m.ss_init[var[tuple(args)]] - m.ss_fin[var[tuple(args)]]))
    var.dev_exp = Expression(*index_sets, rule=dev_exp_rule)

    m.add_component(varname+'_dev'+'_exp', var.dev_exp)

    #for index in var_dev:
    #    print(type(var.dev_exp[index].expr))

    # not sure the expressions (constant terms) are correct, but solve this later
    # next is to debug this, write function to access expression of each variable,
    # and write function to walk expression tree and replace

def update_constraints(fs):
    m = fs.MB_fuel
    for con in m.component_objects(Constraint):
        # get Expression from constraint
        # update expression
        # create new constraint
        # deactivate old constraint

        conname = con.local_name

        if con.dim() == 0:
            index_sets = []
        elif con.dim() == 1:
            index_sets = [con.index_set()]
        elif con.dim() >= 2:
            index_sets = list(con.index_set().set_tuple)

        def con_expr_rule(m, *args):
            idx = tuple(args)
            if idx in con:
                return con[tuple(args)].expr
            else: 
                return None
        con.exp = Expression(*index_sets, rule=con_expr_rule)

        m.add_component(conname+'_exp', con.exp)

        for index in con:
            # walk expression tree of con.exp[index]
            print(con.exp[index].name)
            e = con.exp[index].expr
            replace_variables(e)
            #walk_tree(e)
            # replace nodes
            # deactivate con[index]
        # create new indexed constraint from con.exp

        #con.exp.pprint()

        break
