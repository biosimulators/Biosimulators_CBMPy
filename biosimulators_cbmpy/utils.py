""" Utilities for working with CBMPy

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-01-02
:Copyright: 2021, BioSimulators Team
:License: MIT
"""

from .data_model import SOLVERS, OPTIMIZATION_METHODS
from biosimulators_utils.report.data_model import VariableResults
from biosimulators_utils.sedml.data_model import Variable  # noqa: F401
from biosimulators_utils.utils.core import validate_str_value, parse_value
import numpy
import types  # noqa: F401

__all__ = [
    'apply_algorithm_change_to_simulation_module_method_args',
    'apply_variables_to_simulation_module_method_args',
    'get_simulation_method_args',
    'validate_variables',
    'get_results_of_variables',
    'get_default_solver_module_function_args',
]


def apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args):
    """ Set the value of an argument of a simulation method based on a SED
    algorithm parameter change

    Args:
        method_props (:obj:`dict`): properties of the simulation method
        argument_change (:obj:`AlgorithmParameterChange`): algorithm parameter change
        model (:obj:`cbmpy.CBModel.Model`): model
        module_method_args (:obj:`dict`): dictionary representing the desired simulation function,
            its parent module, and the desired keyword arguments to the function

    Raises:
        :obj:`NotImplementedError`: if the simulation method doesn't support the parameter
        :obj:`ValueError`: if the new value is not a valid value of the parameter
    """
    arg_kisao_id = argument_change.kisao_id

    parameter_props = method_props['parameters'].get(arg_kisao_id, None)
    if parameter_props is None:
        msg = "`{}` is not a parameter of {} ({}). {} suppports the following parameters:\n  - {}".format(
            arg_kisao_id, method_props['name'], method_props['kisao_id'], method_props['name'],
            '\n  - '.join('`{}`: {}'.format(param_kisao_id, method_props['parameters'][param_kisao_id]['name'])
                          for param_kisao_id in sorted(method_props['parameters'].keys()))
        )
        raise NotImplementedError(msg)

    value = argument_change.new_value
    if not validate_str_value(value, parameter_props['type']):
        msg = "`{}` is not a valid value for parameter {} ({}) of {} ({})".format(
            value, parameter_props['name'], arg_kisao_id,
            method_props['name'], method_props['kisao_id'])
        raise ValueError(msg)

    parsed_value = parse_value(value, parameter_props['type'])

    if arg_kisao_id == 'KISAO_0000553':
        solver_name = argument_change.new_value.upper()
        if solver_name not in parameter_props['enum']:
            msg = "`{}` is not a supported solver for {} ({}). The following solvers (KISAO_0000553) are available for {}:\n  - {}".format(
                argument_change.new_value, method_props['name'], method_props['kisao_id'], method_props['name'],
                '\n  - '.join('`' + solver + '`' for solver in sorted(parameter_props['enum'])))
            raise NotImplementedError(msg)

        module_method_args['solver'] = SOLVERS[solver_name]
        if not SOLVERS[solver_name]['module']:
            raise ModuleNotFoundError('{} solver ({}) is not available.'.format(argument_change.new_value, arg_kisao_id))

    elif arg_kisao_id == 'KISAO_0000552':
        if argument_change.new_value.lower() not in OPTIMIZATION_METHODS:
            msg = ("`{}` is not a supported optimization method. "
                   "The following optimization methods (KISAO_0000552) are available:\n  - {}").format(
                argument_change.new_value, '\n  - '.join('`' + name + '`' for name in sorted(OPTIMIZATION_METHODS)))
            raise NotImplementedError(msg)
        module_method_args['optimization_method'] = argument_change.new_value.lower()

    elif arg_kisao_id == 'KISAO_0000534':
        rxn_ids = set(reaction.id for reaction in model.reactions)
        desired_rxn_ids = set(parsed_value)

        invalid_rxn_ids = desired_rxn_ids.difference(rxn_ids)
        if invalid_rxn_ids:
            msg = (
                'Some of the values of {} ({}) of {} ({}) are not SBML ids of reactions:\n  - {}\n\n'
                'The values of {} should be drawn from the following list of the SMBL ids of the reactions of the model:\n  - {}'
            ).format(
                parameter_props['name'], arg_kisao_id,
                method_props['name'], method_props['kisao_id'],
                '\n  - '.join('`' + value + '`' for value in sorted(invalid_rxn_ids)),
                parameter_props['name'],
                '\n  - '.join('`' + rxn_id + '`' for rxn_id in sorted(rxn_ids)),
            )
            raise ValueError(msg)

        parsed_value = sorted(desired_rxn_ids)
        module_method_args['args'][parameter_props['arg_name']] = parsed_value

    elif arg_kisao_id == 'KISAO_0000531':
        if parsed_value < 0. or parsed_value > 1.:
            msg = 'Optimum fraction (KISAO_0000531) must be greater than or equal to 0 and less than or equal to 1, not {}.'.format(
                parsed_value)
            raise ValueError(msg)

        parsed_value *= 100.

        module_method_args['args'][parameter_props['arg_name']] = parsed_value


def apply_variables_to_simulation_module_method_args(target_sbml_id_map, method_props, variables, solver_method_args):
    """ Encode the desired output variables into arguments to simulation methods

    Args:
        target_sbml_id_map (:obj:`dict` of :obj:`str` to :obj:`str`): dictionary that maps each XPath to the
            SBML id of the corresponding model object
        method_props (:obj:`dict`): properties of the simulation method
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        solver_method_args (:obj:`dict`): keyword arguments for the simulation method
    """
    if method_props['function_suffix'] == 'FluxVariabilityAnalysis':
        selected_reactions = set()
        for variable in variables:
            selected_reactions.add(target_sbml_id_map[variable.target])
        solver_method_args['selected_reactions'] = sorted(selected_reactions)


def get_simulation_method_args(method_props, module_method_args):
    """ Setup the simulation method and its keyword arguments

    Args:
        method_props (:obj:`dict`): properties of the simulation method
        module_method_args (:obj:`dict`): dictionary representing the desired simulation function,
            its parent module, and the desired keyword arguments to the function

    Returns:
        :obj:`tuple`:

            * :obj:`types.FunctionType`: simulation method
            * :obj:`dict`: keyword arguments for the simulation method
    """
    solver_props = module_method_args['solver']
    solver_module = solver_props['module']
    solver_method_name = solver_props['function_prefix'] + '_' + method_props['function_suffix']
    solver_method = getattr(solver_module, solver_method_name)
    if module_method_args['optimization_method']:
        opt_method = module_method_args['optimization_method']
        module_method_args['args']['method'] = solver_props['optimization_methods'].get(opt_method, None)
        if not module_method_args['args']['method']:
            msg = ("`{}` is not a supported optimization method of {}. "
                   "{} supports the following optimization methods (KISAO_0000552):\n  - {}").format(
                opt_method,
                solver_props['name'],
                solver_props['name'],
                '\n  - '.join('`' + name + '`' for name in sorted(solver_props['optimization_methods'].keys())))
            raise NotImplementedError(msg)

    solver_method_args = dict(**module_method_args['args'], **method_props['default_args'])

    return solver_method, solver_method_args


def validate_variables(model, method, variables, target_sbml_id_map, target_sbml_fbc_id_map, sbml_fbc_uri):
    """ Validate the desired output variables of a simulation

    Args:
        model (:obj:`cbmpy.CBModel.Model`): model
        method (:obj:`dict`): properties of desired simulation method
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        target_sbml_id_map (:obj:`dict` of :obj:`str` to :obj:`str`): dictionary that maps each XPath to the
            SBML id of the corresponding model object
        target_sbml_fbc_id_map (:obj:`dict` of :obj:`str` to :obj:`str`): dictionary that maps each XPath to the
            SBML-FBC id of the corresponding model object
        sbml_fbc_uri (:obj:`str`): URI for SBML FBC package
    """
    possible_target_results_path_map = set()
    for variable_pattern in method['variables']:
        for sbml_id, fbc_id, attr, _, _ in variable_pattern['get_target_results_paths'](model):
            possible_target_results_path_map.add((sbml_id, fbc_id, attr))

    invalid_symbols = set()
    invalid_targets = set()
    for variable in variables:
        if variable.symbol:
            invalid_symbols.add(variable.symbol)

        else:
            valid = True

            target = variable.target
            variable_target_id = target_sbml_id_map.get(target, None)
            variable_target_fbc_id = target_sbml_fbc_id_map.get(target, None)
            target_attr = target.partition('/@')[2] or None
            if target_attr:
                target_ns, _, target_attr = target_attr.rpartition(':')
                if target_ns and variable.target_namespaces.get(target_ns, None) != sbml_fbc_uri:
                    valid = False

            if not valid or (variable_target_id, variable_target_fbc_id, target_attr) not in possible_target_results_path_map:
                invalid_targets.add(variable.target)

    if invalid_symbols:
        msg = "{} ({}) doesn't support variables with symbols".format(
            method['name'], method['kisao_id'])
        raise NotImplementedError(msg)

    if invalid_targets:
        msg = (
            "{} ({}) doesn't support variables with the following target XPaths:\n  - {}\n\n"
            "The targets of variables should match one of the following patterns of XPaths:\n  - {}"
        ).format(
            method['name'], method['kisao_id'],
            '\n  - '.join(sorted('`' + target + '`' for target in invalid_targets)),
            '\n  - '.join(sorted('{}: `{}`'.format(
                variable_pattern['description'], variable_pattern['target'])
                for variable_pattern in method['variables']))
        )
        raise ValueError(msg)


def get_results_paths_for_variables(model, method_props, variables, target_sbml_id_map, target_sbml_fbc_id_map):
    """ Get the path to results for the desired variables

    Args:
        model (:obj:`cbmpy.CBModel.Model`): model
        method_props (:obj:`dict`): properties of desired simulation method
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        target_sbml_id_map (:obj:`dict` of :obj:`str` to :obj:`str`): dictionary that maps each XPath to the
            SBML id of the corresponding model object
        target_sbml_fbc_id_map (:obj:`dict` of :obj:`str` to :obj:`str`): dictionary that maps each XPath to the
            SBML-FBC id of the corresponding model object

    Returns:
        :obj:`dict`: path to results of desired variables
    """
    possible_target_results_path_map = {}
    for variable_pattern in method_props['variables']:
        for sbml_id, fbc_id, attr, result_type, result_name in variable_pattern['get_target_results_paths'](model):
            possible_target_results_path_map[(sbml_id, fbc_id, attr)] = (result_type, result_name)

    target_results_path_map = {}
    for variable in variables:
        target = variable.target
        variable_target_id = target_sbml_id_map[target]
        variable_target_fbc_id = target_sbml_fbc_id_map[target]
        target_attr = target.partition('/@')[2].rpartition(':')[2] or None
        target_results_path_map[variable.target] = possible_target_results_path_map[(
            variable_target_id, variable_target_fbc_id, target_attr)]

    return target_results_path_map


def get_results_of_variables(target_results_path_map, method_props, solver,
                             variables, model, solution):
    """ Get the results of the desired variables

    Args:
        target_results_path_map (:obj:`dict`): path to results of desired variables
        method_props (:obj:`dict`): properties of desired simulation method
        solver (:obj:`dict`): dictionary representing the desired simulation function,
            its parent module, and the desired keyword arguments to the function
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        model (:obj:`cbmpy.CBModel.Model`): model
        solution (:obj:`object`): solution of method

    Returns:
        :obj:`VariableResults`: the results of desired variables
    """
    all_values = method_props['get_results'](method_props, solver, model, solution)

    variable_results = VariableResults()
    for variable in variables:
        result_type, result_name = target_results_path_map[variable.target]
        result = all_values[result_type].get(result_name, numpy.nan)
        variable_results[variable.id] = numpy.array(result)

    return variable_results


def get_default_solver_module_function_args(algorithm=None):
    """ Get the default solver and its default arguments for an algorithm

    Args:
        algorithm (:obj:`str`, optional): KiSAO id of algorithm

    Returns:
        :obj:`dict`: default solver and default values of its arguments
    """
    return {
        'solver': SOLVERS['GLPK'] if algorithm != 'KISAO_0000554' else SOLVERS['CPLEX'],
        'optimization_method': None,
        'args': {
            'quiet': True,
        }
    }
