"""Startup regressions for the exact generic-grid failure from launcher.log."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from grid_validation import validate_phase_grid, validate_configuration_grids
from launcher import configuration, PRESETS
from mixture_kinetics import MixtureSlab, ROOT
from particle_kinetics import PhaseGrid


class GridInputValidation(unittest.TestCase):
    def setUp(self):
        self.grid=dict(energy_min_eV=1e-4,energy_max_eV=200.,energy_points=128,angular_points=8)

    def test_all_invalid_values_are_named_without_repair(self):
        g=dict(energy_min_eV=0,energy_max_eV=-1,energy_points=8,angular_points=1);before=deepcopy(g)
        with self.assertRaises(ValueError) as error:validate_phase_grid(g,'heavy_phase_grid')
        for text in ['heavy_phase_grid','energy_min_eV=0','energy_max_eV=-1','energy_points=8','angular_points=1']:
            self.assertIn(text,str(error.exception))
        self.assertEqual(g,before)

    def test_nan_infinity_bool_strings_and_fractional_counts_are_rejected(self):
        for key,values in [('energy_min_eV',[float('nan'),float('inf'),True,'0.001']),('energy_max_eV',[float('nan'),float('inf'),True,'200']),('energy_points',[16.5,'128',False]),('angular_points',[3,4.5,'8'])]:
            for value in values:
                with self.subTest(key=key,value=value),self.assertRaisesRegex(ValueError,key):
                    PhaseGrid(dict(self.grid,**{key:value}))

    def test_numeric_boundary_must_exceed_lower_boundary(self):
        with self.assertRaisesRegex(ValueError,'больше'):validate_phase_grid(dict(self.grid,energy_max_eV=1e-5))

    def test_auto_is_input_only_and_never_changes_the_source(self):
        g=dict(self.grid,energy_max_eV='auto');cfg={'schema':'MixtureKineticSlab-2.1','phase_grid':g};before=deepcopy(cfg)
        validate_configuration_grids(cfg);self.assertEqual(cfg,before)
        with self.assertRaisesRegex(ValueError,'auto'):PhaseGrid(g)

    def test_failure_precedes_atomic_data_loading(self):
        cfg={'schema':'MixtureKineticSlab-2.1','phase_grid':dict(self.grid,angular_points=1)}
        with patch('mixture_kinetics.MixtureData') as data:
            with self.assertRaisesRegex(ValueError,'angular_points=1'):MixtureSlab(cfg)
            data.assert_not_called()

    def test_launcher_names_the_same_invalid_field(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'gas_probe.json';p.write_text(json.dumps({'schema':'MixtureKineticSlab-2.1','phase_grid':dict(self.grid,energy_points=8)}),encoding='utf-8')
            # Invalid legacy inputs must still open in the editor for repair.
            self.assertEqual(configuration(p)[0]['phase_grid']['energy_points'],8)
            with self.assertRaisesRegex(ValueError,'energy_points=8'):configuration(p,validate=True)

    def test_first_physics_preset_is_registered_and_structurally_valid(self):
        name='configs/study_01_decay.json';self.assertIn(name,PRESETS)
        cfg,engine=configuration(ROOT/name);self.assertEqual(engine,'mixture_kinetics.py');self.assertEqual(cfg['phase_grid'],self.grid)
        self.assertEqual(cfg['gas']['components'],[{'species':'H_0:L003','density_cm3':1e8}])
