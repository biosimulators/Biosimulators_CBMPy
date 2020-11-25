""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from Biosimulations_utils.simulation.data_model import SteadyStateSimulation, SimulationResultsFormat  # noqa: F401
from Biosimulations_utils.simulator.utils import exec_simulations_in_archive
import cbmpy
import numpy
import json
import os


__all__ = ['exec_combine_archive', 'exec_simulation']

KISAO_SOLVERS = {
    'CPLEX': {
        'id': 'CPLEX',
        'module': cbmpy.CBCPLEX,
        'function_prefix': 'cplx',
        'methods': {
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
        'id': 'GLPK',
        'module': cbmpy.CBGLPK,
        'function_prefix': 'glpk',
        'methods': {
            'simplex': 's',
            'interior': 'i',
            'exact': 'e',
        },
    },
}

KISAO_ALGORITHMS_PARAMETERS_MAP = {
    'KISAO_0000437': {
        'function_suffix': 'analyzeModel',
        'parameters': {
        },
    },
    'KISAO_0000528': {
        'function_suffix': 'MinimizeSumOfAbsFluxes',
        'parameters': {
            'KISAO_0000534': {
                'arg_name': 'selected_reactions',
                'parser': json.loads,
            },
            'KISAO_0000531': {
                'arg_name': 'optPercentage',
                'parser': float,
            },
        }
    },
    'KISAO_0000554': {
        'function_suffix': 'MinimizeNumActiveFluxes',
        'parameters': {
            'KISAO_0000534': {
                'arg_name': 'selected_reactions',
                'parser': json.loads,
            },
            'KISAO_0000531': {
                'arg_name': 'optPercentage',
                'parser': float,
            },
        }
    },
    'KISAO_0000526': {
        'function_suffix': 'FluxVariabilityAnalysis',
        'parameters': {
            'KISAO_0000534': {
                'arg_name': 'selected_reactions',
                'parser': json.loads,
            },
            'KISAO_0000531': {
                'arg_name': 'optPercentage',
                'parser': float,
            },
        }
    }
}


def exec_combine_archive(archive_file, out_dir):
    """ Execute the SED tasks defined in a COMBINE archive and save the outputs

    Args:
        archive_file (:obj:`str`): path to COMBINE archive
        out_dir (:obj:`str`): directory to store the outputs of the tasks
    """
    exec_simulations_in_archive(archive_file, exec_simulation, out_dir, apply_model_changes=True)


def exec_simulation(model_filename, model_sed_urn, simulation, working_dir, out_filename, out_format):
    ''' Execute a simulation and save its results

    Args:
       model_filename (:obj:`str`): path to the model
       model_sed_urn (:obj:`str`): SED URN for the format of the model (e.g., `urn:sedml:language:sbml`)
       simulation (:obj:`SteadyStateSimulation`): simulation
       working_dir (:obj:`str`): directory of the SED-ML file
       out_filename (:obj:`str`): path to save the results of the simulation
       out_format (:obj:`SimulationResultsFormat`): format to save the results of the simulation (e.g., `HDF5`)
    '''
    # check that model is encoded in SBML
    if model_sed_urn != "urn:sedml:language:sbml":
        raise NotImplementedError("Model language with URN '{}' is not supported".format(model_sed_urn))

    # check that simulation is a time course simulation
    if not isinstance(simulation, SteadyStateSimulation):
        raise NotImplementedError('{} is not supported'.format(simulation.__class__.__name__))

    # check that model parameter changes have already been applied (because handled by :obj:`exec_simulations_in_archive`)
    if simulation.model_parameter_changes:
        raise NotImplementedError('Model parameter changes are not supported')

    # check that the desired output format is supported
    if out_format != SimulationResultsFormat.HDF5:
        raise NotImplementedError("Simulation results format '{}' is not supported".format(out_format))

    # Read the model located at `os.path.join(working_dir, model_filename)` in the format
    # with the SED URN `model_sed_urn`.
    model = cbmpy.CBRead.readSBML3FBC(model_filename)

    # Load the algorithm specified by `simulation.algorithm` and parameter values
    algorithm = KISAO_ALGORITHMS_PARAMETERS_MAP.get(simulation.algorithm.kisao_term.id, None)
    if algorithm is None:
        raise NotImplementedError(
            "Algorithm with KiSAO id '{}' is not supported".format(simulation.algorithm.kisao_term.id))

    solver = KISAO_SOLVERS.get('glpk')
    solver_method = None
    solver_args = {
        'with_reduced_costs': True,
        'return_lp_obj': True,
        'quiet': True,
    }
    for parameter_change in simulation.algorithm_parameter_changes:
        if parameter_change.parameter.kisao_term.id == 'KISAO_0000553':
            solver = KISAO_SOLVERS.get(parameter_change.value.upper(), None)
            if solver is None:
                raise NotImplementedError(
                    "Solver with name '{}' is not supported".format(parameter_change.value))
        elif parameter_change.parameter.kisao_term.id == 'KISAO_0000552':
            solver_method = parameter_change.value
        else:
            parameter = algorithm['parameters'].get(parameter_change.parameter.kisao_term.id, None)
            if parameter is None:
                raise NotImplementedError(
                    "Parameter '{}' is not supported for algorithm '{}'".format(
                        parameter_change.parameter.kisao_term.id, simulation.algorithm.kisao_term.id))
            try:
                solver_args[parameter['arg_name']] = parameter['parser'](parameter_change.value)
            except:
                raise ValueError("'{}' is not a valid value of parameter '{}' of '{}'".format(
                    parameter_change.value, parameter_change.parameter.kisao_term.id, simulation.algorithm.kisao_term.id))

    solver_function = getattr(solver['module'], solver['function_prefix'] + '_' + algorithm['function_suffix'])
    if solver_method:
        solver_method_arg = solver['methods'].get(solver_method.lower(), None)
        if solver_method_arg is None:
            raise NotImplementedError(
                "Solver method with name '{}' is not supported".format(solver_method))
        solver_args['method'] = solver_method_arg

    # Simulate the model from `simulation.start_time` to `simulation.end_time`
    solution = solver_function(model, **solver_args)

    # throw error if status isn't optimal
    status = getattr(solver['module'], solver['function_prefix'] + '_' + 'getSolutionStatus')(solution)
    if status != 'LPS_OPT':
        raise ValueError("A solution could not be found. The solver status was '{}'.".format(
            status))

    # get results
    obj_id = model.getActiveObjective().id
    obj_value = model.getObjFuncValue()

    rxn_fluxes = model.getReactionValues()
    rxn_reduced_costs = getattr(solver['module'], solver['function_prefix'] + '_' + 'getReducedCosts')(solution)

    if solver['function_prefix'] != 'cplx':
        raise NotImplementedError("'{}' solver does not support the calculation of shadow prices".format(solver['id']))
    species_shadow_prices = getattr(solver['module'], solver['function_prefix'] + '_' + 'getShadowPrices')(solution)

    # Save a report of the results of the simulation with `simulation.num_time_points` time points
    # beginning at `simulation.output_start_time` to `out_filename` in `out_format` format.
    # This should save all of the variables specified by `simulation.model.variables`.

    pass
