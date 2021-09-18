""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import KISAO_ALGORITHMS_PARAMETERS_MAP
from .utils import (apply_algorithm_change_to_simulation_module_method_args,
                    apply_variables_to_simulation_module_method_args,
                    get_simulation_method_args, validate_variables,
                    get_results_paths_for_variables, get_results_of_variables,
                    get_default_solver_module_function_args)
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.config import get_config, Config  # noqa: F401
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog, StandardOutputErrorCapturerLevel  # noqa: F401
from biosimulators_utils.model_lang.sbml.utils import get_package_namespace as get_sbml_package_namespace
from biosimulators_utils.viz.data_model import VizFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults, SedDocumentResults  # noqa: F401
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, ModelAttributeChange, SteadyStateSimulation,  # noqa: F401
                                                  Variable)
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.exec import exec_sed_doc as base_exec_sed_doc
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.utils.core import raise_errors_warnings
from biosimulators_utils.warnings import warn, BioSimulatorsWarning
from biosimulators_utils.xml.utils import get_namespaces_for_xml_doc
from kisao.data_model import AlgorithmSubstitutionPolicy, ALGORITHM_SUBSTITUTION_POLICY_LEVELS
from kisao.utils import get_preferred_substitute_algorithm_by_ids
from lxml import etree
import cbmpy
import copy

__all__ = [
    'exec_sedml_docs_in_combine_archive',
    'exec_sed_doc',
    'exec_sed_task',
    'preprocess_sed_task',
]


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir, config=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            * :obj:`SedDocumentResults`: results
            * :obj:`CombineArchiveLog`: log
    """
    return exec_sedml_docs_in_archive(exec_sed_doc, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      config=config)


def exec_sed_doc(doc, working_dir, base_out_path, rel_out_path=None,
                 apply_xml_model_changes=True,
                 log=None, indent=0, pretty_print_modified_xml_models=False,
                 log_level=StandardOutputErrorCapturerLevel.c, config=None):
    """ Execute the tasks specified in a SED document and generate the specified outputs

    Args:
        doc (:obj:`SedDocument` or :obj:`str`): SED document or a path to SED-ML file which defines a SED document
        working_dir (:obj:`str`): working directory of the SED document (path relative to which models are located)

        base_out_path (:obj:`str`): path to store the outputs

            * CSV: directory in which to save outputs to files
              ``{base_out_path}/{rel_out_path}/{report.id}.csv``
            * HDF5: directory in which to save a single HDF5 file (``{base_out_path}/reports.h5``),
              with reports at keys ``{rel_out_path}/{report.id}`` within the HDF5 file

        rel_out_path (:obj:`str`, optional): path relative to :obj:`base_out_path` to store the outputs
        apply_xml_model_changes (:obj:`bool`, optional): if :obj:`True`, apply any model changes specified in the SED-ML file before
            calling :obj:`task_executer`.
        log (:obj:`SedDocumentLog`, optional): log of the document
        indent (:obj:`int`, optional): degree to indent status messages
        pretty_print_modified_xml_models (:obj:`bool`, optional): if :obj:`True`, pretty print modified XML models
        log_level (:obj:`StandardOutputErrorCapturerLevel`, optional): level at which to log output
        config (:obj:`Config`, optional): BioSimulators common configuration
        simulator_config (:obj:`SimulatorConfig`, optional): tellurium configuration

    Returns:
        :obj:`tuple`:

            * :obj:`ReportResults`: results of each report
            * :obj:`SedDocumentLog`: log of the document
    """
    return base_exec_sed_doc(exec_sed_task, doc, working_dir, base_out_path,
                             rel_out_path=rel_out_path,
                             apply_xml_model_changes=apply_xml_model_changes,
                             log=log,
                             indent=indent,
                             pretty_print_modified_xml_models=pretty_print_modified_xml_models,
                             log_level=log_level,
                             config=config)


def exec_sed_task(task, variables, preprocessed_task=None, log=None, config=None):
    ''' Execute a task and save its results

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        preprocessed_task (:obj:`dict`, optional): preprocessed information about the task, including possible
            model changes and variables. This can be used to avoid repeatedly executing the same initialization
            for repeated calls to this method.
        log (:obj:`TaskLog`, optional): log for the task
        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            :obj:`VariableResults`: results of variables
            :obj:`TaskLog`: log

    Raises:
        :obj:`ValueError`: if the task or an aspect of the task is not valid, or the requested output variables
            could not be recorded
        :obj:`NotImplementedError`: if the task is not of a supported type or involves an unsuported feature
    '''
    config = config or get_config()

    if config.LOG and not log:
        log = TaskLog()

    if preprocessed_task is None:
        preprocessed_task = preprocess_sed_task(task, variables, config=config)

    model = task.model

    # Read the model
    cbmpy_model = preprocessed_task['model']['model']

    # modify model
    if model.changes:
        raise_errors_warnings(validation.validate_model_change_types(model.changes, (ModelAttributeChange,)),
                              error_summary='Changes for model `{}` are not supported.'.format(model.id))
        model_change_setter_map = preprocessed_task['model']['model_change_setter_map']
        for change in model.changes:
            new_value = float(change.new_value)
            model_change_setter_map[change.target](new_value)

    # Set up simulation function and its keyword arguments
    variable_target_sbml_id_map = preprocessed_task['model']['variable_target_sbml_id_map']
    method_props = preprocessed_task['simulation']['method_props']
    simulation_method_args = copy.copy(preprocessed_task['simulation']['method_args'])
    apply_variables_to_simulation_module_method_args(variable_target_sbml_id_map, method_props, variables, simulation_method_args)

    # Simulate the model
    simulation_method = preprocessed_task['simulation']['method']
    solution = simulation_method(cbmpy_model, **simulation_method_args)

    # throw error if status isn't optimal
    module_method_args = preprocessed_task['simulation']['module_method_args']
    method_props['raise_if_simulation_error'](module_method_args, solution)

    # get the result of each variable
    variable_results = get_results_of_variables(
        preprocessed_task['model']['variable_target_results_path_map'],
        method_props, module_method_args['solver'],
        variables, cbmpy_model, solution)

    # log action
    if config.LOG:
        log.algorithm = preprocessed_task['simulation']['algorithm_kisao_id']
        log.simulator_details = {
            'method': simulation_method.__module__ + '.' + simulation_method.__name__,
            'arguments': simulation_method_args,
        }

    # return the result of each variable and log
    return variable_results, log


def preprocess_sed_task(task, variables, config=None):
    """ Preprocess a SED task, including its possible model changes and variables. This is useful for avoiding
    repeatedly initializing tasks on repeated calls of :obj:`exec_sed_task`.

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`dict`: preprocessed information about the task
    """
    config = config or get_config()

    model = task.model
    sim = task.simulation

    if config.VALIDATE_SEDML:
        raise_errors_warnings(validation.validate_task(task),
                              error_summary='Task `{}` is invalid.'.format(task.id))
        raise_errors_warnings(validation.validate_model_language(model.language, ModelLanguage.SBML),
                              error_summary='Language for model `{}` is not supported.'.format(model.id))
        raise_errors_warnings(validation.validate_model_change_types(model.changes, (ModelAttributeChange,)),
                              error_summary='Changes for model `{}` are not supported.'.format(model.id))
        raise_errors_warnings(*validation.validate_model_changes(model),
                              error_summary='Changes for model `{}` are invalid.'.format(model.id))
        raise_errors_warnings(validation.validate_simulation_type(sim, (SteadyStateSimulation, )),
                              error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
        raise_errors_warnings(*validation.validate_simulation(sim),
                              error_summary='Simulation `{}` is invalid.'.format(sim.id))
        raise_errors_warnings(*validation.validate_data_generator_variables(variables),
                              error_summary='Data generator variables for task `{}` are invalid.'.format(task.id))

    model_etree = etree.parse(model.source)
    model_change_sbml_id_map = validation.validate_target_xpaths(
        model.changes, model_etree, attr='id')
    variable_target_sbml_id_map = validation.validate_target_xpaths(
        variables, model_etree, attr='id')
    namespaces = get_namespaces_for_xml_doc(model_etree)
    sbml_fbc_prefix, sbml_fbc_uri = get_sbml_package_namespace('fbc', namespaces)
    variable_target_sbml_fbc_id_map = validation.validate_target_xpaths(
        variables,
        model_etree,
        attr={
            'namespace': {
                'prefix': sbml_fbc_prefix,
                'uri': sbml_fbc_uri,
            },
            'name': 'id',
        }
    )

    # Read the model
    cbmpy_model = cbmpy.CBRead.readSBML3FBC(model.source)

    # preprocess model changes
    model_change_setter_map = {}
    sbml_id_reaction_map = {rxn.id: rxn for rxn in cbmpy_model.reactions}
    invalid_changes = []
    for change in model.changes:
        target = change.target
        sbml_id = model_change_sbml_id_map[target]
        reaction = sbml_id_reaction_map.get(sbml_id, None)
        setter = None
        if reaction:
            attr = target.partition('/@')[2]
            ns, _, attr = attr.partition(':')
            ns = change.target_namespaces.get(ns, None)

            if ns == sbml_fbc_uri:
                if attr == 'lowerFluxBound':
                    setter = reaction.setLowerBound
                elif attr == 'upperFluxBound':
                    setter = reaction.setUpperBound
        if setter:
            model_change_setter_map[change.target] = setter
        else:
            invalid_changes.append(change.target)

    if invalid_changes:
        valid_changes = []
        for reaction in cbmpy_model.reactions:
            valid_changes.append(
                "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='{}']/@fbc:lowerFluxBound".format(reaction.id))
            valid_changes.append(
                "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='{}']/@fbc:upperFluxBound".format(reaction.id))

        msg = 'The following changes are invalid:\n  {}\n\nThe following targets are valid:\n  {}'.format(
            '\n  '.join(sorted(invalid_changes)),
            '\n  '.join(sorted(valid_changes)),
        )
        raise ValueError(msg)

    # Set up the algorithm specified by :obj:`task.simulation.algorithm.kisao_id`
    alg_substitution_policy = get_algorithm_substitution_policy(config=config)
    exec_kisao_id = get_preferred_substitute_algorithm_by_ids(
        sim.algorithm.kisao_id, KISAO_ALGORITHMS_PARAMETERS_MAP.keys(),
        substitution_policy=alg_substitution_policy)
    method_props = KISAO_ALGORITHMS_PARAMETERS_MAP[exec_kisao_id]

    # Set up the the parameters of the algorithm
    module_method_args = get_default_solver_module_function_args(exec_kisao_id)
    if exec_kisao_id == sim.algorithm.kisao_id:
        for change in sim.algorithm.changes:
            try:
                apply_algorithm_change_to_simulation_module_method_args(method_props, change, cbmpy_model, module_method_args)
            except NotImplementedError as exception:
                if (
                    ALGORITHM_SUBSTITUTION_POLICY_LEVELS[alg_substitution_policy]
                    > ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.NONE]
                ):
                    warn('Unsuported algorithm parameter `{}` was ignored:\n  {}'.format(
                        change.kisao_id, str(exception).replace('\n', '\n  ')),
                        BioSimulatorsWarning)
                else:
                    raise
            except ValueError as exception:
                if (
                    ALGORITHM_SUBSTITUTION_POLICY_LEVELS[alg_substitution_policy]
                    > ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.NONE]
                ):
                    warn('Unsuported value `{}` for algorithm parameter `{}` was ignored:\n  {}'.format(
                        change.new_value, change.kisao_id, str(exception).replace('\n', '\n  ')),
                        BioSimulatorsWarning)
                else:
                    raise

    # validate selected solver is available
    if not module_method_args['solver']['module']:
        raise ModuleNotFoundError('{} solver is not available.'.format(module_method_args['solver']['name']))

    # Setup simulation function and its keyword arguments
    simulation_method, simulation_method_args = get_simulation_method_args(method_props, module_method_args)

    # validate and preprocess variables
    validate_variables(cbmpy_model, method_props, variables, variable_target_sbml_id_map, variable_target_sbml_fbc_id_map, sbml_fbc_uri)
    variable_target_results_path_map = get_results_paths_for_variables(
        cbmpy_model, method_props, variables, variable_target_sbml_id_map, variable_target_sbml_fbc_id_map)

    # return preprocessed information
    return {
        'model': {
            'model': cbmpy_model,
            'model_change_setter_map': model_change_setter_map,
            'variable_target_sbml_id_map': variable_target_sbml_id_map,
            'variable_target_results_path_map': variable_target_results_path_map,
        },
        'simulation': {
            'algorithm_kisao_id': exec_kisao_id,
            'method': simulation_method,
            'method_props': method_props,
            'method_args': simulation_method_args,
            'module_method_args': module_method_args,
        }
    }
