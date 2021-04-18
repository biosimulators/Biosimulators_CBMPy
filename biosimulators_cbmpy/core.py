""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import KISAO_ALGORITHMS_PARAMETERS_MAP, DEFAULT_SOLVER_MODULE_FUNCTION_ARGS
from .utils import (apply_algorithm_change_to_simulation_module_method_args,
                    apply_variables_to_simulation_module_method_args,
                    get_simulation_method_args, validate_variables,
                    get_results_of_variables)
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog  # noqa: F401
from biosimulators_utils.plot.data_model import PlotFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults  # noqa: F401
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, SteadyStateSimulation,  # noqa: F401
                                                  Variable)
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.exec import exec_sed_doc
from biosimulators_utils.utils.core import raise_errors_warnings
from biosimulators_utils.xml.utils import get_namespaces_for_xml_doc
from lxml import etree
import cbmpy
import copy
import functools

__all__ = [
    'exec_sedml_docs_in_combine_archive',
    'exec_sed_task',
]


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                       report_formats=None, plot_formats=None,
                                       bundle_outputs=None, keep_individual_outputs=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        report_formats (:obj:`list` of :obj:`ReportFormat`, optional): report format (e.g., csv or h5)
        plot_formats (:obj:`list` of :obj:`PlotFormat`, optional): report format (e.g., pdf)
        bundle_outputs (:obj:`bool`, optional): if :obj:`True`, bundle outputs into archives for reports and plots
        keep_individual_outputs (:obj:`bool`, optional): if :obj:`True`, keep individual output files

    Returns:
        :obj:`CombineArchiveLog`: log
    """
    sed_doc_executer = functools.partial(exec_sed_doc, exec_sed_task)
    return exec_sedml_docs_in_archive(sed_doc_executer, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      report_formats=report_formats,
                                      plot_formats=plot_formats,
                                      bundle_outputs=bundle_outputs,
                                      keep_individual_outputs=keep_individual_outputs)


def exec_sed_task(task, variables, log=None):
    ''' Execute a task and save its results

    Args:
       task (:obj:`Task`): task
       variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
       log (:obj:`TaskLog`, optional): log for the task

    Returns:
        :obj:`tuple`:

            :obj:`VariableResults`: results of variables
            :obj:`TaskLog`: log

    Raises:
        :obj:`ValueError`: if the task or an aspect of the task is not valid, or the requested output variables
            could not be recorded
        :obj:`NotImplementedError`: if the task is not of a supported type or involves an unsuported feature
    '''
    log = log or TaskLog()

    model = task.model
    sim = task.simulation

    raise_errors_warnings(validation.validate_task(task),
                          error_summary='Task `{}` is invalid.'.format(task.id))
    raise_errors_warnings(validation.validate_model_language(task.model.language, ModelLanguage.SBML),
                          error_summary='Language for model `{}` is not supported.'.format(model.id))
    raise_errors_warnings(validation.validate_model_change_types(task.model.changes, ()),
                          error_summary='Changes for model `{}` are not supported.'.format(model.id))
    raise_errors_warnings(validation.validate_model_changes(task.model),
                          error_summary='Changes for model `{}` are invalid.'.format(model.id))
    raise_errors_warnings(validation.validate_simulation_type(task.simulation, (SteadyStateSimulation, )),
                          error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
    raise_errors_warnings(validation.validate_simulation(task.simulation),
                          error_summary='Simulation `{}` is invalid.'.format(sim.id))
    raise_errors_warnings(validation.validate_data_generator_variables(variables),
                          error_summary='Data generator variables for task `{}` are invalid.'.format(task.id))
    target_x_paths_ids = validation.validate_variable_xpaths(
        variables, task.model.source, attr='id')
    namespaces = get_namespaces_for_xml_doc(etree.parse(task.model.source))
    target_x_paths_fbc_ids = validation.validate_variable_xpaths(
        variables,
        task.model.source,
        attr={
            'namespace': {
                'prefix': 'fbc',
                'uri': namespaces['fbc'],
            },
            'name': 'id',
        }
    )

    # Read the model
    model = cbmpy.CBRead.readSBML3FBC(task.model.source)

    # Set up the algorithm specified by :obj:`task.simulation.algorithm.kisao_id`
    simulation = task.simulation
    algorithm_kisao_id = simulation.algorithm.kisao_id
    method_props = KISAO_ALGORITHMS_PARAMETERS_MAP.get(algorithm_kisao_id, None)
    if method_props is None:
        msg = "".join([
            "Algorithm with KiSAO id `{}` is not supported. ".format(algorithm_kisao_id),
            "Algorithm must have one of the following KiSAO ids:\n  - {}".format('\n  - '.join(
                '`{}`: {}'.format(kisao_id, method_props['name'])
                for kisao_id, method_props in KISAO_ALGORITHMS_PARAMETERS_MAP.items())),
        ])
        raise NotImplementedError(msg)

    # Set up the the parameters of the algorithm
    module_method_args = copy.copy(DEFAULT_SOLVER_MODULE_FUNCTION_ARGS)
    module_method_args['args'] = copy.copy(module_method_args['args'])
    for change in simulation.algorithm.changes:
        apply_algorithm_change_to_simulation_module_method_args(method_props, change, model, module_method_args)

    # validate variables
    validate_variables(method_props, variables)

    # set keyword arguments based on desired outputs
    apply_variables_to_simulation_module_method_args(target_x_paths_ids, method_props, variables, module_method_args)

    # Setup simulation function and its keyword arguments
    simulation_method, simulation_method_args = get_simulation_method_args(method_props, module_method_args)

    # Simulate the model
    solution = simulation_method(model, **simulation_method_args)

    # throw error if status isn't optimal
    method_props['raise_if_simulation_error'](module_method_args, solution)

    # get the result of each variable
    variable_results = get_results_of_variables(target_x_paths_ids, target_x_paths_fbc_ids,
                                                method_props, module_method_args['solver'],
                                                variables, model, solution)

    # log action
    log.algorithm = algorithm_kisao_id
    log.simulator_details = {
        'method': simulation_method.__module__ + '.' + simulation_method.__name__,
        'arguments': simulation_method_args,
    }

    # return the result of each variable and log
    return variable_results, log
