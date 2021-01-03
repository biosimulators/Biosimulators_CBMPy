from biosimulators_cbmpy.data_model import SOLVERS, KISAO_ALGORITHMS_PARAMETERS_MAP, DEFAULT_SOLVER_MODULE_FUNCTION_ARGS
from biosimulators_cbmpy.utils import (apply_algorithm_change_to_simulation_module_method_args,
                                       get_simulation_method_kw_args,
                                       validate_variables, get_results_of_variables)
from biosimulators_utils.sedml.data_model import AlgorithmParameterChange, DataGeneratorVariable
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

    def test_apply_algorithm_change_to_simulation_module_method_args(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000526']
        model = mock.Mock()
        module_method_args = copy.copy(DEFAULT_SOLVER_MODULE_FUNCTION_ARGS)
        module_method_args['kw_args'] = copy.copy(module_method_args['kw_args'])
        expected_module_method_args = copy.copy(DEFAULT_SOLVER_MODULE_FUNCTION_ARGS)
        expected_module_method_args['kw_args'] = copy.copy(expected_module_method_args['kw_args'])

        # other parameters
        argument_change = AlgorithmParameterChange(
            kisao_id='KISAO_0000531',
            new_value='0.99',
        )
        apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)
        expected_module_method_args['kw_args']['optPercentage'] = 0.99 * 100
        self.assertEqual(module_method_args, expected_module_method_args)

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
        with self.assertRaisesRegex(ValueError, 'not available'):
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
        argument_change = AlgorithmParameterChange(
            kisao_id='KISAO_0000534',
            new_value='["A", "B"]',
        )
        apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)
        expected_module_method_args['kw_args']['selected_reactions'] = ['A', 'B']
        self.assertEqual(module_method_args, expected_module_method_args)

        argument_change.new_value = '["A", "B", "C"]'
        with self.assertRaisesRegex(ValueError, 'not SBML ids of reactions'):
            apply_algorithm_change_to_simulation_module_method_args(method_props, argument_change, model, module_method_args)

    def test_get_simulation_method_kw_args(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000437']

        module_method_args = copy.copy(DEFAULT_SOLVER_MODULE_FUNCTION_ARGS)
        module_method_args['kw_args'] = copy.copy(module_method_args['kw_args'])

        simulation_method, simulation_method_kw_args = get_simulation_method_kw_args(method_props, module_method_args)

        self.assertEqual(simulation_method, cbmpy.CBGLPK.glpk_analyzeModel)
        self.assertEqual(simulation_method_kw_args, {
            'quiet': True,
            'with_reduced_costs': True,
            'return_lp_obj': True,
        })

        # set optimization method
        module_method_args['optimization_method'] = 'simplex'
        simulation_method, simulation_method_kw_args = get_simulation_method_kw_args(method_props, module_method_args)
        self.assertEqual(simulation_method, cbmpy.CBGLPK.glpk_analyzeModel)
        self.assertEqual(simulation_method_kw_args, {
            'method': 's',
            'quiet': True,
            'with_reduced_costs': True,
            'return_lp_obj': True,
        })

        # set optimization method
        module_method_args['optimization_method'] = 'auto'
        with self.assertRaisesRegex(NotImplementedError, 'not a supported optimization method'):
            simulation_method, simulation_method_kw_args = get_simulation_method_kw_args(method_props, module_method_args)

    def test_validate_variables(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000437']
        variables = [
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']/@value"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:type='maximize']/@value"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective/@value"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@flux"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@reducedCost"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@metaid='R_ACALD']/@flux"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction/@flux"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']/@shadowPrice"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@metaid='M_13dpg_c']/@shadowPrice"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species/@shadowPrice"),
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species"),
        ]
        validate_variables(method_props, variables)

        variables = [
            DataGeneratorVariable(symbol='urn:sedml:symbol:time'),
        ]
        with self.assertRaises(NotImplementedError):
            validate_variables(method_props, variables)

        variables = [
            DataGeneratorVariable(target="/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment[@id='c']")
        ]
        with self.assertRaises(ValueError):
            validate_variables(method_props, variables)

    def test_get_results_of_variables(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000437']
        solver = SOLVERS['GLPK']
        variables = [
            DataGeneratorVariable(id='obj',
                                  target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='obj']/@value"),
            DataGeneratorVariable(id='inactive_obj',
                                  target="/sbml:sbml/sbml:model/fbc:listOfObjectives/fbc:objective[@fbc:id='inactive_obj']/@value"),
            DataGeneratorVariable(id='R_ACALD_flux',
                                  target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@flux"),
            DataGeneratorVariable(id='R_ACALD_reduced_cost',
                                  target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_ACALD']/@reducedCost"),
            DataGeneratorVariable(id='M_13dpg_c_shadow_price',
                                  target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='M_13dpg_c']/@shadowPrice"),
            DataGeneratorVariable(id='M_2pg_c_shadow_price',
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
            species=[mock.Mock(id='M_13dpg_c'), mock.Mock(id='M_2pg_c')],
            reactions=[mock.Mock(id='R_ACALD')],
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
        result = get_results_of_variables(target_to_id, target_to_fbc_id, method_props, solver, variables, model, solution)
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
        result = get_results_of_variables(target_to_id, target_to_fbc_id, method_props, solver, variables, model, solution)
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
        result = get_results_of_variables(target_to_id, target_to_fbc_id, method_props, solver, variables, model, solution)
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
            DataGeneratorVariable(id='R_PGI_min_flux',
                                  target="/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction[@id='R_PGI']/@minFlux"),
            DataGeneratorVariable(id='R_PGI_max_flux',
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
        result = get_results_of_variables(target_to_id, target_to_fbc_id, method_props, solver, variables, model, solution)
        numpy.testing.assert_allclose(result['R_PGI_min_flux'], numpy.array(-15.))
        numpy.testing.assert_allclose(result['R_PGI_max_flux'], numpy.array(25.))

    def test_raise_if_simulation_error(self):
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000437']
        module_method_args = copy.copy(DEFAULT_SOLVER_MODULE_FUNCTION_ARGS)
        module_method_args['kw_args'] = copy.copy(module_method_args['kw_args'])

        # optimal FBA solution
        solution = mock.Mock(status='opt')
        method_props['raise_if_simulation_error'](module_method_args, solution)

        # infeasible FBA solution
        solution.status = 'infeas'
        with self.assertRaisesRegex(ValueError, ''):
            method_props['raise_if_simulation_error'](module_method_args, solution)

        # FVA
        method_props = KISAO_ALGORITHMS_PARAMETERS_MAP['KISAO_0000526']
        method_props['raise_if_simulation_error'](None, None)
