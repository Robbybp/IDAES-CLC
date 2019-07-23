from pyomo.environ import (Var, SimpleVar, Expression, Constraint, Suffix,
        Param, value)
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
        #print('node:', node)
        if (node.__class__ in native_numeric_types or 
                isinstance(node, NumericConstant) or 
                isinstance(node, Param)):
            # why would node ever be None?
            return True, node

        try:
            if isinstance(node.parent_component(), Param):
                return True, node,
        except AttributeError:
            pass

        if node.is_variable_type():
            # probably includes expressions including variables as well
            # in which case should return false...
            if isinstance(node, ExpressionBase):
                print('node is an instance of ExpressionBase')
            try: 
                # okay for now, but node might not have parent component but still have dev exp
                # ^ this doesn't seem to make sense, but I have no proof
                # SimpleVar.parent_component() will return the SimpleVar
                # (Thanks John!)
                if node.parent_component().has_dev_exp == True:
                    #print('exp:', node.parent_component().dev_exp[node.index()])
                    #print('replacing')
                    return True, node.parent_component().dev_exp[node.index()].expr
                else:
                    return True, node
            except AttributeError:
                # catch for variables with no deviation expression
                return True, node
                
        return False, node

class TestWalker(SimpleExpressionVisitor):

    def __init__(self):
        super(SimpleExpressionVisitor, self).__init__()

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
#    pdb.set_trace()
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

def create_scale_values(var, fs, fs_ss_init, fs_ss_fin, tol=1e-8):
    m = fs.MB_fuel
    m_ss_init = fs_ss_init.MB_fuel
    m_ss_fin = fs_ss_fin.MB_fuel
    time = m.t
    varname = var.local_name

    if not is_indexed_by(var, time):
        var.has_dev_exp = False
        if isinstance(Var, SimpleVar):
            if var.value != 0:
                m.scaling_factor[var] = var.value
        else:
            for index in var:
                if var[index].value != 0:
                    m.scaling_factor[var[index]] = var[index].value
        return

    #if isinstance(var, DerivativeVar):
    # not obvious what to do about derivatives...
    # scale, same as variable, or not?
    # really z-derivatives should be scaled same as any other variable
    #    var.has_dev_exp = False
    #    return
    # ^will assume that the function is not called for time derivatives

    var_ss_init = eval('m_ss_init.'+varname)
    var_ss_fin = eval('m_ss_fin.'+varname)

    if var.index_set().dimen == 1:
        index_sets = [var.index_set()]
    elif var.index_set().dimen >= 2:
        index_sets = list(var.index_set().set_tuple)

    # need to remember which var data objects have deviation variables
    #var.dev = Var(*index_sets,initialize=

    var.has_dev_exp = False
    var_fixed = True
    var_effectively_fixed = True

    for index in var:
        # this loop assigns suffix values for each variable
        nt_index = get_non_time_index(var, index, time)
        m.ss_init[var[index]] = var_ss_init[nt_index].value
        m.ss_fin[var[index]] = var_ss_fin[nt_index].value

        if var[index].value == None:
            continue

        if var[index].fixed == True:
            pass
        elif var[index].fixed == False:
            var_fixed = False

        if var == m.dG_fluxdz:
            print(index)
            #pdb.set_trace()

        if abs(var_ss_init[nt_index].value - 
               var_ss_fin[nt_index].value) < tol:
            m.scaling_factor[var[index]] = var_ss_init[nt_index].value

        elif abs(var_ss_init[nt_index].value - 
                 var_ss_fin[nt_index].value) >= tol:
            m.scaling_factor[var[index]] = (var_ss_init[nt_index].value -
                                            var_ss_fin[nt_index].value)
            var.has_dev_exp = True
            var_effectively_fixed = False
            
    if var_fixed == True or var_effectively_fixed == True:
        var.has_dev_exp = False
        # but for variables that are fixed, still need (want) a scaling factor
        # ^ provided above
        return

    def dev_init_rule(m, *args):
        # can a variable indexed by a single set accept a tuple as a subscript?
        # I hope so...
        if not var[tuple(args)].value == None:
            return var[tuple(args)].value - m.ss_fin[var[tuple(args)]]
    # will be initialized to zero for "non-dynamic" variables
    # will this cause any problems? should this be initialized to 
    # the var's value?
    var.dev = Var(*index_sets, initialize=dev_init_rule)

    m.add_component(varname + '_dev', var.dev)
    var_dev = eval('m.'+varname+'_dev')

    for index in var_dev:
        if var_dev[index].value != 0 and var_dev[index].value != None:
            m.scaling_factor[var_dev[index]] = m.scaling_factor[var[index]]

    def dev_exp_rule(m, *args):
        if var[tuple(args)].value == None:
            return
        if abs(m.ss_init[var[tuple(args)]] - 
                m.ss_fin[var[tuple(args)]]) < tol:
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

def update_constraint(con, fs):
    # to test:
    # eq_p2 for a constraint that will not have Cg replaced
    # eq_c6 or eq_c4 for constraits that will have Cg replaced
    m = fs.MB_fuel
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
        # component indexed by a single set can be subscripted by a tuple,
        # but tuple(idx) in component will return False
        idx = None
        if len(args) == 1:
            idx = args[0]
        elif len(args) > 1:
            idx = tuple(args)
        if idx is None:
            return con.expr
        if idx in con:
            return con[tuple(args)].expr
        else: 
            return None
    con.exp = Expression(*index_sets, rule=con_expr_rule)
    # when a continuous set is present, does this behave as expected?

    m.add_component(conname+'_exp', con.exp)

    idx_exp_map = {}
    # ^ entire purpose of following loop is to construct this map
    for index in con:
        #if con == m.eq_a1:
        #    print(index, con[index].expr.to_string())
        # walk expression tree of con.exp[index]
        #print(con.exp[index].name)
        e = con.exp[index].expr
        # need to access this expression from con[index] somehow
        idx_exp_map[index] = replace_variables(e)
        # ^ what happens to e here? -TODO investigate-

        #print(idx_exp_map[index])
        #print(value(idx_exp_map[index]))
        # deactivate con[index]
    # create new indexed constraint from con.exp

    # this seems to give the desired result
    #print(idx_exp_map[(0.00062,0)])

    def dev_con_rule(m, *args):
        idx = None
        if len(args) == 1:
            idx = args[0]
        elif len(args) > 1:
            idx = tuple(args)
        if idx in con:
            #print(idx_exp_map[idx])
            # ^ this does not seem to give the desired result
            # had wrong variable as the subscript
            return idx_exp_map[idx]
        else:
            return Constraint.Skip
    con.dev_con = Constraint(*index_sets, rule=dev_con_rule)
    m.add_component(conname+'_dev', con.dev_con)

    # EqualityExpression (ExpressionBase) has no attribute 'pprint()'
    # but IndexedExpression does
    #con.exp.pprint()
    #con.dev_exp.pprint()
