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
import re
import types  # noqa: F401

__all__ = [
    'apply_algorithm_change_to_simulation_module_method_args',
    'apply_variables_to_simulation_module_method_args',
    'get_simulation_method_args',
    'validate_variables',
    'get_results_of_variables',
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
            raise ValueError('{} solver ({}) is not available.'.format(argument_change.new_value, arg_kisao_id))

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
        parsed_value *= 100.

        module_method_args['args'][parameter_props['arg_name']] = parsed_value


def apply_variables_to_simulation_module_method_args(target_x_paths_ids, method_props, variables, module_method_args):
    """ Encode the desired output variables into arguments to simulation methods

    Args:
        target_x_paths_ids (:obj:`dict` of :obj:`str` to :obj:`str`): dictionary that maps each XPath to the
            SBML id of the corresponding model object
        method_props (:obj:`dict`): properties of the simulation method
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        module_method_args (:obj:`dict`): dictionary representing the desired simulation function,
            its parent module, and the desired keyword arguments to the function
    """
    if method_props['kisao_id'] == 'KISAO_0000526':
        selected_reactions = set()
        for variable in variables:
            selected_reactions.add(target_x_paths_ids[variable.target])
        module_method_args['args']['selected_reactions'] = sorted(selected_reactions)


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


def validate_variables(method, variables):
    """ Validate the desired output variables of a simulation

    Args:
        method (:obj:`dict`): properties of desired simulation method
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
    """
    invalid_symbols = set()
    invalid_targets = set()
    for variable in variables:
        if variable.symbol:
            invalid_symbols.add(variable.symbol)

        else:
            valid = False
            for variable_pattern in method['variables']:
                if re.match(variable_pattern['target'], variable.target):
                    valid = True
                    break

            if not valid:
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


def get_results_of_variables(target_x_paths_ids, target_x_paths_fbc_ids, method_props, solver,
                             variables, model, solution):
    """ Get the results of the desired variables

    Args:
        target_x_paths_ids (:obj:`dict` of :obj:`str` to :obj:`str`): dictionary that maps each XPath to the
            SBML id of the corresponding model object
        target_x_paths_fbc_ids (:obj:`dict` of :obj:`str` to :obj:`str`): dictionary that maps each XPath to the
            SBML-FBC id of the corresponding model object
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
        target = variable.target
        for variable_pattern in method_props['variables']:
            if re.match(variable_pattern['target'], target):
                variable_target_id = target_x_paths_ids[target]
                variable_target_fbc_id = target_x_paths_fbc_ids[target]
                result = variable_pattern['get_result'](variable_target_id,
                                                        variable_target_fbc_id, all_values)

                break

        variable_results[variable.id] = numpy.array(result)

    return variable_results
