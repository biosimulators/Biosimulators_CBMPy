from biosimulators_cbmpy.data_model import SOLVERS, KISAO_ALGORITHMS_PARAMETERS_MAP
from biosimulators_cbmpy.utils import (apply_algorithm_change_to_simulation_module_method_args,
                                       apply_variables_to_simulation_module_method_args,
                                       get_simulation_method_args,
                                       validate_variables, get_results_paths_for_variables, get_results_of_variables,
                                       get_default_solver_module_function_args)
from biosimulators_utils.sedml.data_model import AlgorithmParameterChange, Variable
from numpy import nan
from unittest import mock
import cbmpy
import copy
import numpy
import numpy.testing
import os
import unittest


class UtilsTestCase(unittest.TestCase):
    MODEL_FILENAME = os.path.join(os.path.dirname(__file__), 'fixtures', 'textbook.xml')
    NAMESPACES = {
        'sbml': 'http://www.sbml.org/sbml/level3/version1/core',
        'fbc': 'http://www.sbml.org/sbml/level3/version1/fbc/version2',
    }

    def test_apply_algorithm_change_to_simulation_module_method_args(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000526']
        model = mock.Mock()
        module_method_args = get_default_solver_module_function_args()
        expected_module_method_args = get_default_solver_module_function_args()

        # other parameters
        argument_change = AlgorithmParameterChange(
            kisao_id='KISAO_0000531',
            new_value='0.99',
        )
        apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)
        expected_module_method_args['args']['optPercentage'] = 0.99 * 100
        self.assertEqual(module_method_args, expected_module_method_args)

        argument_change.new_value = 2
        with self.assertRaisesRegex(ValueError, 'less than or equal'):
            apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)

        argument_change.new_value = 'text'
        with self.assertRaisesRegex(ValueError, 'not a valid value'):
            apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)

        argument_change.kisao_id = 'KISAO_0000001'
        with self.assertRaisesRegex(NotImplementedError, 'not a parameter of'):
            apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)

        # solver
        SOLVERS['CPLEX']['module'] = True
        argument_change = AlgorithmParameterChange(
            kisao_id='KISAO_0000553',
            new_value='CPLEX',
        )
        apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)
        expected_module_method_args['solver'] = SOLVERS['CPLEX']
        self.assertEqual(module_method_args, expected_module_method_args)

        argument_change.new_value = 'unsupported'
        with self.assertRaisesRegex(NotImplementedError, 'not a supported solver for'):
            apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)

        SOLVERS['CPLEX']['module'] = None
        argument_change.new_value = 'CPLEX'
        with self.assertRaisesRegex(ModuleNotFoundError, 'not available'):
            apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)
        SOLVERS['CPLEX']['module'] = True

        # optimization method
        argument_change = AlgorithmParameterChange(
            kisao_id='KISAO_0000552',
            new_value='PRIMAL',
        )
        apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)
        expected_module_method_args['optimization_method'] = 'primal'
        self.assertEqual(module_method_args, expected_module_method_args)

        argument_change.new_value = 'unsupported'
        with self.assertRaisesRegex(NotImplementedError, 'not a supported optimization method'):
            apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)

        # reaction list
        model.reactions = [
            mock.Mock(id='A'),
            mock.Mock(id='B'),
        ]
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000528']
        argument_change = AlgorithmParameterChange(
            kisao_id='KISAO_0000534',
            new_value='["A", "B"]',
        )
        apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)
        expected_module_method_args['args']['selected_reactions'] = ['A', 'B']
        self.assertEqual(module_method_args, expected_module_method_args)

        argument_change.new_value = '["A", "B", "C"]'
        with self.assertRaisesRegex(ValueError, 'not SBML ids of reactions'):
            apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)

    def test_apply_variables_to_simulation_module_method_args(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000526']
        variables = [
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='A']/@minFlux"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='A']/@maxFlux"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='B']/@minFlux"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='C']/@maxFlux"),
        ]
        target_x_paths_ids = {
            variables[0].target: 'A',
            variables[1].target: 'A',
            variables[2].target: 'B',
            variables[3].target: 'C',
        }

        # FVA
        module_method_args = {'args': {}}
        expected_module_method_args = {'args': {'selected_reactions': ['A', 'B', 'C']}}
        apply_variables_to_simulation_module_method_args(target_x_paths_ids, method_props, variables, module_method_args['args'])
        self.assertEqual(module_method_args['args'], expected_module_method_args['args'])

        # FBA
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000437']
        module_method_args = {'args': {}}
        expected_module_method_args = {'args': {}}
        apply_variables_to_simulation_module_method_args(target_x_paths_ids, method_props, variables, module_method_args['args'])
        self.assertEqual(module_method_args['args'], expected_module_method_args['args'])

    def test_get_simulation_method_args(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000437']

        module_method_args = get_default_solver_module_function_args()

        simulation_method, simulation_method_args = get_simulation_method_args(method_props, module_method_args)

        self.assertEqual(simulation_method, cbmpy.CBGLPK.glpk_analyzeModel)
        self.assertEqual(simulation_method_args, {
            'quiet': True,
            'with_reduced_costs': True,
            'return_lp_obj': True,
        })

        # set optimization method
        module_method_args['optimization_method'] = 'simplex'
        simulation_method, simulation_method_args = get_simulation_method_args(method_props, module_method_args)
        self.assertEqual(simulation_method, cbmpy.CBGLPK.glpk_analyzeModel)
        self.assertEqual(simulation_method_args, {
            'method': 's',
            'quiet': True,
            'with_reduced_costs': True,
            'return_lp_obj': True,
        })

        # set optimization method
        module_method_args['optimization_method'] = 'auto'
        with self.assertRaisesRegex(NotImplementedError, 'not a supported optimization method'):
            simulation_method, simulation_method_args = get_simulation_method_args(method_props, module_method_args)

    def test_validate_variables(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000437']

        model = cbmpy.CBRead.readSBML3FBC(self.MODEL_FILENAME)
        variable_target_sbml_id_map = {
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']/@value": None,
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:type='maximize']/@value": None,
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']": None,
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective/@value": None,
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective": None,
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@flux": 'R_ACALD',
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@reducedCost": 'R_ACALD',
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@metaid='R_ACALD']/@flux": 'R_ACALD',
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']": 'R_ACALD',
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction/@flux": None,
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction": None,
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']/@shadowPrice": 'M_13dpg_c',
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@metaid='M_13dpg_c']/@shadowPrice": 'M_13dpg_c',
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']": 'M_13dpg_c',
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species/@shadowPrice": None,
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species": None,
        }
        variable_target_sbml_fbc_id_map = {
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']/@value": 'obj',
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:type='maximize']/@value": 'obj',
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']": 'obj',
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective/@value": None,
            "/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective": None,
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@flux": None,
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@reducedCost": None,
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@metaid='R_ACALD']/@flux": None,
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']": None,
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction/@flux": None,
            "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction": None,
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']/@shadowPrice": None,
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@metaid='M_13dpg_c']/@shadowPrice": None,
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']": None,
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species/@shadowPrice": None,
            "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species": None,
        }      
        sbml_fbc_uri = self.NAMESPACES['fbc']

        variables = [
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']/@value"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:type='maximize']/@value"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@flux"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@reducedCost"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@metaid='R_ACALD']/@flux"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']/@shadowPrice"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@metaid='M_13dpg_c']/@shadowPrice"),
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']"),
        ]
        validate_variables(model, method_props, variables, variable_target_sbml_id_map, variable_target_sbml_fbc_id_map, sbml_fbc_uri)

        variables = [
            Variable(symbol='urn:sedml:symbol:time'),
        ]
        with self.assertRaises(NotImplementedError):
            validate_variables(model, method_props, variables, variable_target_sbml_id_map, variable_target_sbml_fbc_id_map, sbml_fbc_uri)

        variables = [
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment[@id='c']")
        ]
        with self.assertRaises(ValueError):
            validate_variables(model, method_props, variables, variable_target_sbml_id_map, variable_target_sbml_fbc_id_map, sbml_fbc_uri)

        variables = [
            Variable(target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment[@id='c']/@sbml:name")
        ]
        with self.assertRaises(ValueError):
            validate_variables(model, method_props, variables, variable_target_sbml_id_map, variable_target_sbml_fbc_id_map, sbml_fbc_uri)

    def test_get_results_of_variables(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000437']
        solver = SOLVERS['GLPK']
        variables = [
            Variable(id='obj',
                     target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']/@value"),
            Variable(id='inactive_obj',
                     target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='inactive_obj']/@value"),
            Variable(id='R_ACALD_flux',
                     target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@flux"),
            Variable(id='R_ACALD_reduced_cost',
                     target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@reducedCost"),
            Variable(id='M_13dpg_c_shadow_price',
                     target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']/@shadowPrice"),
            Variable(id='M_2pg_c_shadow_price',
                     target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_2pg_c']/@shadowPrice"),

        ]

        target_to_id = {
            variables[0].target: None,
            variables[1].target: None,
            variables[2].target: 'R_ACALD',
            variables[3].target: 'R_ACALD',
            variables[4].target: 'M_13dpg_c',
            variables[5].target: 'M_2pg_c',
        }
        target_to_fbc_id = {
            variables[0].target: 'obj',
            variables[1].target: 'inactive_obj',
            variables[2].target: None,
            variables[3].target: None,
            variables[4].target: None,
            variables[5].target: None,
        }

        model = cbmpy.CBRead.readSBML3FBC(self.MODEL_FILENAME)
        model = mock.Mock(
            getActiveObjective=lambda: mock.Mock(id='obj'),
            getObjFuncValue=lambda: 0.8739215069684909,
            getReactionValues=lambda: {'R_ACALD': 1.250555e-12},
            objectives=[mock.Mock(id='obj'), mock.Mock(id='inactive_obj')],
            species=[mock.Mock(id='M_13dpg_c'), mock.Mock(id='M_2pg_c')],
            reactions=[mock.Mock(id='R_ACALD'), mock.Mock(id='R_PGI'), mock.Mock(id='R_PPC')],
        )

        # FBA, GLPK
        # solution = cbmpy.CBGLPK.glpk_analyzeModel(model, with_reduced_costs=True, return_lp_obj=True, quiet=True)
        solver = {
            'name': 'GLPK',
            'function_prefix': 'glpk',
            'module': mock.Mock(
                glpk_getReducedCosts=lambda solution: {'R_ACALD': 0.},
            )
        }
        solution = None
        target_results_path_map = get_results_paths_for_variables(model, method_props, variables, target_to_id, target_to_fbc_id)
        result = get_results_of_variables(target_results_path_map, method_props, solver, variables, model, solution)
        self.assertEqual(set(result.keys()), set(var.id for var in variables))
        numpy.testing.assert_allclose(result['obj'], numpy.array(0.8739215069684909))
        numpy.testing.assert_allclose(result['inactive_obj'], numpy.array(nan))
        numpy.testing.assert_allclose(result['R_ACALD_flux'], numpy.array(1.250555e-12), rtol=1e-6, atol=1e-8)
        numpy.testing.assert_allclose(result['R_ACALD_reduced_cost'], numpy.array(0.,), rtol=1e-6, atol=1e-8)
        numpy.testing.assert_allclose(result['M_13dpg_c_shadow_price'], numpy.array(nan))
        numpy.testing.assert_allclose(result['M_2pg_c_shadow_price'], numpy.array(nan))

        # pFBA, CPLEX
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000528']
        # solver = SOLVERS['CPLEX']
        # solution = cbmpy.CBCPLEX.cplx_analyzeModel(model, with_reduced_costs=True, return_lp_obj=True, quiet=True)
        solver = {
            'name': 'CPLEX',
            'function_prefix': 'cplx',
            'module': mock.Mock(
                cplx_getReducedCosts=lambda solution: {'R_ACALD': 0.},
                cplx_getShadowPrices=lambda solution: {
                    'M_13dpg_c': (-49.330418, 0, 1.0),
                    'M_2pg_c': (-35., 0, 65.),
                },
            )
        }
        solution = None
        target_results_path_map = get_results_paths_for_variables(model, method_props, variables, target_to_id, target_to_fbc_id)
        result = get_results_of_variables(target_results_path_map, method_props, solver, variables, model, solution)
        self.assertEqual(set(result.keys()), set(var.id for var in variables))
        numpy.testing.assert_allclose(result['obj'], numpy.array(0.8739215069684909))
        numpy.testing.assert_allclose(result['inactive_obj'], numpy.array(nan))
        numpy.testing.assert_allclose(result['R_ACALD_flux'], numpy.array(1.250555e-12), rtol=1e-6, atol=1e-8)
        numpy.testing.assert_allclose(result['R_ACALD_reduced_cost'], numpy.array(0.,), rtol=1e-6, atol=1e-8)
        numpy.testing.assert_allclose(result['M_13dpg_c_shadow_price'], numpy.array(-49.330418), rtol=1e-6, atol=1e-8)
        numpy.testing.assert_allclose(result['M_2pg_c_shadow_price'], numpy.array(65.), rtol=1e-6, atol=1e-8)

        # pFBA, CPLEX
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000554']
        # solver = SOLVERS['CPLEX']
        # solution = cbmpy.CBCPLEX.cplx_analyzeModel(model, with_reduced_costs=True, return_lp_obj=True, quiet=True)
        solver = {
            'name': 'CPLEX',
            'function_prefix': 'cplx',
            'module': mock.Mock()
        }
        solution = (None, None)
        target_results_path_map = get_results_paths_for_variables(model, method_props, variables, target_to_id, target_to_fbc_id)
        result = get_results_of_variables(target_results_path_map, method_props, solver, variables, model, solution)
        self.assertEqual(set(result.keys()), set(var.id for var in variables))
        numpy.testing.assert_allclose(result['obj'], numpy.array(0.8739215069684909))
        numpy.testing.assert_allclose(result['inactive_obj'], numpy.array(nan))
        numpy.testing.assert_allclose(result['R_ACALD_flux'], numpy.array(1.250555e-12), rtol=1e-6, atol=1e-8)
        numpy.testing.assert_allclose(result['R_ACALD_reduced_cost'], numpy.array(nan))
        numpy.testing.assert_allclose(result['M_13dpg_c_shadow_price'], numpy.array(nan))
        numpy.testing.assert_allclose(result['M_2pg_c_shadow_price'], numpy.array(nan))

        # FVA
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000526']
        solver = SOLVERS['GLPK']
        variables = [
            Variable(id='R_PGI_min_flux',
                     target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_PGI']/@minFlux"),
            Variable(id='R_PGI_max_flux',
                     target_namespaces=self.NAMESPACES,
                     target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_PGI']/@maxFlux"),
        ]
        target_to_id = {
            variables[0].target: 'R_PGI',
            variables[1].target: 'R_PGI',
        }
        target_to_fbc_id = {
            variables[0].target: None,
            variables[1].target: None,
        }
        solution = (
            numpy.array(
                [
                    [nan] * 7,
                    [nan] * 2 + [-15, 25.] + [nan] * 3,
                    [nan] * 7,
                ],
            ),
            ['R_ACALD', 'R_PGI', 'R_PPC']
        )
        target_results_path_map = get_results_paths_for_variables(model, method_props, variables, target_to_id, target_to_fbc_id)
        result = get_results_of_variables(target_results_path_map, method_props, solver, variables, model, solution)
        numpy.testing.assert_allclose(result['R_PGI_min_flux'], numpy.array(-15.))
        numpy.testing.assert_allclose(result['R_PGI_max_flux'], numpy.array(25.))

    def test_raise_if_simulation_error(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000437']
        module_method_args = get_default_solver_module_function_args()

        # optimal FBA solution
        model = cbmpy.CBRead.readSBML3FBC(self.MODEL_FILENAME)
        solution = cbmpy.CBGLPK.glpk_analyzeModel(model, return_lp_obj=True)
        method_props['raise_if_simulation_error'](module_method_args, solution)

        # infeasible FBA solution
        model = cbmpy.CBRead.readSBML3FBC(self.MODEL_FILENAME)
        for bound in model.flux_bounds:
            bound.setValue(1000)
        solution = cbmpy.CBGLPK.glpk_analyzeModel(model, return_lp_obj=True)
        with self.assertRaisesRegex(ValueError, ''):
            method_props['raise_if_simulation_error'](module_method_args, solution)

        # FVA
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000526']
        method_props['raise_if_simulation_error'](None, None)

        model = cbmpy.CBRead.readSBML3FBC(self.MODEL_FILENAME)
        solution = cbmpy.CBGLPK.glpk_FluxVariabilityAnalysis(model)
        method_props['raise_if_simulation_error'](module_method_args, solution)
