# Handles solvers
from __future__ import division
from bunch import Bunch
from six import iteritems

def get_solvers():
    """Returns a dictionary (actually a Bunch) of solver classes

    Returns:
        Bunch: solver classes
    """
    cplex_lib = "/opt/ibm/ILOG/"\
        "CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so"

    solvers = Bunch()
    solvers.local_NLP = SolverClass()
    solvers.local_NLP.solver_list = ['ipopt', 'conopt', 'baron', 'knitro', 'minos', 'snopt']

    # solvers.local_NLP.ipopt.linear_solver = 'ma57'
    # solvers.local_NLP.ipopt.linear_system_scaling = 'mc19'
    # solvers.local_NLP.ipopt.mu_init = '1E-10'
    # solvers.local_NLP.ipopt.bound_push = '1E-10'

    solvers.petsc = Bunch()
    solvers.petsc['-snes_monitor'] = ''
    solvers.petsc['-snes_max_it'] = '20000'
    solvers.petsc['-snes_max_funcs'] = '60000'
    solvers.petsc['-snes_atol'] = '1e-10'
    solvers.petsc['-snes_rtol'] = '1e-8'

    solvers.global_NLP = SolverClass()
    solvers.global_NLP.solver_list = ['baron', 'couenne']
    solvers.global_NLP.baron.CplexLibName = cplex_lib
    solvers.global_NLP.baron.allowipopt = 0
    solvers.global_NLP.baron.EpsA = 0
    solvers.global_NLP.baron.EpsR = 0.001
    solvers.global_NLP.baron.MaxTime = 60

    solvers.global_MINLP = SolverClass()
    solvers.global_MINLP.solver_list = ['baron', 'couenne']
    solvers.global_MINLP.baron.CplexLibName = cplex_lib
    solvers.global_MINLP.baron.allowipopt = 0
    solvers.global_MINLP.baron.EpsA = 0
    solvers.global_MINLP.baron.EpsR = 0.001
    solvers.global_MINLP.baron.MaxTime = 600

    solvers.self_proj_cut_gen = SolverClass()
    solvers.self_proj_cut_gen.solver_list = ['baron', 'couenne']
    solvers.self_proj_cut_gen.baron.CplexLibName = cplex_lib
    solvers.self_proj_cut_gen.baron.allowipopt = 0
    solvers.self_proj_cut_gen.baron.EpsA = 0
    solvers.self_proj_cut_gen.baron.EpsR = 0.001
    solvers.self_proj_cut_gen.baron.MaxTime = 60

    solvers.local_MIP = SolverClass()
    solvers.local_MIP.solver_list = ['gurobi', 'cplex', 'cbc']

    solvers.bounding = SolverClass()
    solvers.bounding.solver_list = ['gurobi', 'cplex', 'cbc']

    solvers.global_MIP = SolverClass()
    solvers.global_MIP.solver_list = ['gurobi', 'cplex', 'cbc']

    return solvers

class SolverError(EnvironmentError):
    pass

class SolverClass(Bunch):
    def __init__(self, *args, **kwargs):
        super(SolverClass, self).__init__(*args, **kwargs)
        self._solver_list = []
        self._reserved = set(['_reserved', '_solver_list', 'solve'])

    def solver_options(self, s_name=''):
        """Generator for the solver options

        Args:
            s_name (str, optional): solver name

        Returns:
            tuple: (option name, option value)
        """
        for key, value in iteritems(self):
            if key not in self._reserved and key not in self.solver_list:
                yield key, value
        if s_name:
            for key, value in iteritems(self[s_name]):
                if key not in self._reserved:
                    yield key, value

    @property
    def solver_list(self):
        return self._solver_list

    @solver_list.setter
    def solver_list(self, value):
        for old_solver in set(self._solver_list) - set(value):
            del self.old_solver
        for new_solver in set(value) - set(self._solver_list):
            self[new_solver] = Bunch()
        self._solver_list = value

    def promote(self, sol):
        """Promote a solver to the start of the solver list

        Args:
            sol (str): solver name

        Returns:
            None
        """
        self.solver_list.remove(sol)
        self.solver_list.insert(0, sol)
