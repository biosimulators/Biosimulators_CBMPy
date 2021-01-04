""" Data model for mapping from KiSAO terms for algorithms and their parameters
to solvers methods and their arguments.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-01-02
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_utils.data_model import ValueType
from numpy import nan
import cbmpy  # noqa: F401
import numpy
try:
    from cbmpy import CBGLPK
except ModuleNotFoundError:  # pragma: no cover  # GLPK should always be installed
    CBGLPK = None
try:
    from cbmpy import CBCPLEX
except ModuleNotFoundError:
    CBCPLEX = None

__all__ = [
    'SOLVERS',
    'OPTIMIZATION_METHODS',
    'FBA_DEPENDENT_VARIABLE_TARGETS',
    'FVA_DEPENDENT_VARIABLE_TARGETS',
    'KISAO_ALGORITHMS_PARAMETERS_MAP',
    'DEFAULT_SOLVER_MODULE_FUNCTION_ARGS',
]

SOLVERS = {
    'CPLEX': {
        'name': 'CPLEX',
        'module': CBCPLEX,
        'function_prefix': 'cplx',
        'optimization_methods': {
            'auto': 'o',
            'primal': 'p',
            'dual': 'd',
            'barrier without crossover': 'b',
            'barrier': 'h',
            'sifting': 's',
            'concurrent': 'c',
        },
    },
    'GLPK': {
        'name': 'GLPK',
        'module': CBGLPK,
        'function_prefix': 'glpk',
        'optimization_methods': {
            'simplex': 's',
            'interior': 'i',
            'exact': 'e',
        },
    },
}

OPTIMIZATION_METHODS = set([
    'auto',
    'primal',
    'dual',
    'barrier without crossover',
    'barrier',
    'sifting',
    'concurrent',
    'simplex',
    'interior',
    'exact',
])

FBA_DEPENDENT_VARIABLE_TARGETS = [
    {
        'description': 'objective value',
        'target_type': 'objective',
        'target': r'^/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective(\[.*?\])?(/@value)?$',
        'get_result':
            lambda el_id, el_fbc_id, solution:
                solution['objective_value'].get(el_fbc_id, nan),
    },
    {
        'description': 'reaction flux',
        'target_type': 'reaction',
        'target': r'^/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction(\[.*?\])?(/@flux)?$',
        'get_result':
            lambda el_id, el_fbc_id, solution:
                solution['reaction_flux'].get(el_id),
    },
    {
        'description': 'reaction reduced cost',
        'target_type': 'reaction',
        'target': r'^/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction(\[.*?\])?/@reducedCost$',
        'get_result':
            lambda el_id, el_fbc_id, solution:
                solution['reaction_reduced_cost'][el_id],
    },
    {
        'description': 'species shadow price',
        'target_type': 'species',
        'target': r'^/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species(\[.*?\])?(/@shadowPrice)?$',
        'get_result':
            lambda el_id, el_fbc_id, solution:
                solution['species_shadow_price'][el_id],
    },
]

FVA_DEPENDENT_VARIABLE_TARGETS = [
    {
        'description': 'minimum reaction flux',
        'target_type': 'reaction',
        'target': r'^/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction(\[.*?\])?/@minFlux?$',
        'get_result':
            lambda el_id, el_fbc_id, solution:
                solution['reaction_min_flux'][el_id],
    },
    {
        'description': 'maximum reaction flux',
        'target_type': 'reaction',
        'target': r'^/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction(\[.*?\])?/@maxFlux?$',
        'get_result':
            lambda el_id, el_fbc_id, solution:
                solution['reaction_max_flux'][el_id],
    },
]


def raise_if_fba_simulation_error(module_method_args, solution):
    """ Raise an error if the solution of an FBA or pFBA simulation is not optimal

    Args:
        module_method_args (:obj:`dict`): dictionary representing the desired simulation function,
            its parent module, and the desired keyword arguments to the function
        solution (:obj:`glpk.LPX` or :obj:`cplex.Cplex`): solution of method

    Raises:
        :obj:`ValueError`: error
    """
    solver_props = module_method_args['solver']
    status_method = getattr(solver_props['module'], solver_props['function_prefix'] + '_getSolutionStatus')
    status = status_method(solution)
    if status not in ['LPS_OPT', 'MILP_OPT']:
        raise ValueError("A solution could not be found. The solver status was `{}`.".format(
            status))


def get_fba_results(method_props, solver, model, solution):
    """ Get the results of an FBA simulation

    Args:
        method_props (:obj:`dict`): properties of desired simulation method
        solver (:obj:`dict`): dictionary representing the desired simulation function,
            its parent module, and the desired keyword arguments to the function
        model (:obj:`cbmpy.CBModel.Model`): model
        solution (:obj:`glpk.LPX` or :obj:`cplex.Cplex`): solution of method

    Returns:
        :obj:`dict`: results of an FBA simulation
    """
    results = {}

    results['objective_value'] = {
        model.getActiveObjective().id: model.getObjFuncValue()
    }

    results['reaction_flux'] = model.getReactionValues()

    if method_props['kisao_id'] in ['KISAO_0000437', 'KISAO_0000528']:
        results['reaction_reduced_cost'] = getattr(solver['module'], solver['function_prefix'] + '_' + 'getReducedCosts')(solution)
    else:
        results['reaction_reduced_cost'] = {reaction.id: numpy.nan for reaction in model.reactions}

    if method_props['kisao_id'] in ['KISAO_0000437', 'KISAO_0000528'] and solver['name'] == 'CPLEX':
        results['species_shadow_price'] = {}
        for key, (lb, rhs, ub) in getattr(solver['module'], solver['function_prefix'] + '_' + 'getShadowPrices')(solution).items():
            if ub > -lb:
                results['species_shadow_price'][key] = ub
            else:
                results['species_shadow_price'][key] = lb

    else:
        results['species_shadow_price'] = {species.id: numpy.nan for species in model.species}

    return results


def get_fva_results(method_props, solver, model, solution):
    """ Get the results of an FVA simulation

    Args:
        method_props (:obj:`dict`): properties of desired simulation method
        solver (:obj:`dict`): dictionary representing the desired simulation function,
            its parent module, and the desired keyword arguments to the function
        model (:obj:`cbmpy.CBModel.Model`): model
        solution (:obj:`tuple`): solution of method

            * :obj:`numpy.ndarray`: 2-D array of FVA results.

               * Rows: reactions.
               * Columns:
                    * ``Reaction`
                    * ``Reduced Costs``
                    * ``Variability Min``
                    * ``Variability Max``
                    * ``abs(Max-Min)``
                    * ``MinStatus``
                    * ``MaxStatus``

            * :obj:`list` of :obj:`str`: reaction id for each row

    Returns:
        :obj:`dict`: results of an FVA simulation
    """
    columns = ['Reaction', 'Reduced Costs', 'Variability Min', 'Variability Max', 'abs(Max-Min)', 'MinStatus', 'MaxStatus']
    values, rxn_ids = solution
    assert values.shape[1] == len(columns)

    i_rxn_min_fluxes = columns.index('Variability Min')
    i_rxn_max_fluxes = columns.index('Variability Max')

    results = {
        'reaction_min_flux': {},
        'reaction_max_flux': {},
    }

    for i_reaction, rxn_id in enumerate(rxn_ids):
        results['reaction_min_flux'][rxn_id] = values[i_reaction, i_rxn_min_fluxes]
        results['reaction_max_flux'][rxn_id] = values[i_reaction, i_rxn_max_fluxes]

    return results


KISAO_ALGORITHMS_PARAMETERS_MAP = {
    'KISAO_0000437': {
        'kisao_id': 'KISAO_0000437',
        'name': 'FBA: flux balance analysis',
        'function_suffix': 'analyzeModel',
        'parameters': {
            'KISAO_0000553': {
                'name': 'solver',
                'type': ValueType.string,
                'enum': ['CPLEX', 'GLPK'],
            },
            'KISAO_0000552': {
                'name': 'optimization method',
                'type': ValueType.string,
                'enum': OPTIMIZATION_METHODS,
            },
        },
        'default_args': {
            'with_reduced_costs': True,
            'return_lp_obj': True,
        },
        'variables': FBA_DEPENDENT_VARIABLE_TARGETS,
        'raise_if_simulation_error': raise_if_fba_simulation_error,
        'get_results': get_fba_results,
    },
    'KISAO_0000528': {
        'kisao_id': 'KISAO_0000528',
        'name': 'pFBA: parsimonious flux balance analysis (minimum sum of absolute fluxes)',
        'function_suffix': 'MinimizeSumOfAbsFluxes',
        'parameters': {
            'KISAO_0000534': {
                'name': 'selected reactions',
                'arg_name': 'selected_reactions',
                'type': ValueType.list,
            },
            'KISAO_0000531': {
                'name': 'optimum percentage',
                'arg_name': 'optPercentage',
                'type': ValueType.float,
            },
            'KISAO_0000553': {
                'name': 'solver',
                'type': ValueType.string,
                'enum': ['CPLEX', 'GLPK'],
            },
            'KISAO_0000552': {
                'name': 'optimization method',
                'type': ValueType.string,
                'enum': OPTIMIZATION_METHODS,
            },
        },
        'default_args': {
            'with_reduced_costs': True,
            'return_lp_obj': True,
        },
        'variables': FBA_DEPENDENT_VARIABLE_TARGETS,
        'raise_if_simulation_error': raise_if_fba_simulation_error,
        'get_results': get_fba_results,
    },
    'KISAO_0000554': {
        'kisao_id': 'KISAO_0000554',
        'name': 'pFBA: parsimonious flux balance analysis (minimum number of active fluxes)',
        'function_suffix': 'MinimizeNumActiveFluxes',
        'parameters': {
            'KISAO_0000534': {
                'name': 'selected reactions',
                'arg_name': 'selected_reactions',
                'type': ValueType.list,
            },
            'KISAO_0000531': {
                'name': 'optimum percentage',
                'arg_name': 'optPercentage',
                'type': ValueType.float,
            },
            'KISAO_0000553': {
                'name': 'solver',
                'type': ValueType.string,
                'enum': ['CPLEX'],
            }
        },
        'default_args': {
            'return_lp_obj': True,
        },
        'variables': FBA_DEPENDENT_VARIABLE_TARGETS,
        'raise_if_simulation_error':
            lambda module_method_args, opt_solution:
                raise_if_fba_simulation_error(module_method_args, opt_solution[1]),
        'get_results':
            lambda method_props, solver, model, opt_solution:
                get_fba_results(method_props, solver, model, opt_solution[1]),
    },
    'KISAO_0000526': {
        'kisao_id': 'KISAO_0000526',
        'name': 'FVA: flux variability analysis',
        'function_suffix': 'FluxVariabilityAnalysis',
        'parameters': {
            'KISAO_0000531': {
                'name': 'optimum percentage',
                'arg_name': 'optPercentage',
                'type': ValueType.float,
            },
            'KISAO_0000553': {
                'name': 'solver',
                'type': ValueType.string,
                'enum': ['CPLEX', 'GLPK'],
            },
            'KISAO_0000552': {
                'name': 'optimization method',
                'type': ValueType.string,
                'enum': OPTIMIZATION_METHODS,
            },
        },
        'default_args': {
        },
        'variables': FVA_DEPENDENT_VARIABLE_TARGETS,
        'raise_if_simulation_error': lambda module_method_args, solution: None,
        'get_results': get_fva_results,
    }
}


DEFAULT_SOLVER_MODULE_FUNCTION_ARGS = {
    'solver': SOLVERS['GLPK'],
    'optimization_method': None,
    'args': {
        'quiet': True,
    }
}
