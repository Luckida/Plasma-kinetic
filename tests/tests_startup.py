"""Bootstrap checks require only the standard library, not numerical packages."""
import json
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

import launch_probe


class StartupChecks(unittest.TestCase):
    def probe(self, version=(3, 11, 0), bits=64, implementation='CPython', tk_error=None):
        with patch.object(launch_probe.sys, 'version_info', version), \
             patch.object(launch_probe.struct, 'calcsize', return_value=bits // 8), \
             patch.object(launch_probe.platform, 'python_implementation', return_value=implementation), \
             patch.object(launch_probe, 'check_tk', **({'side_effect': tk_error} if tk_error else {'return_value': '8.6.13'})):
            return launch_probe.inspect_runtime()

    def test_python_311_64_with_tk_is_accepted(self):
        report = self.probe()
        self.assertTrue(report['ok']);self.assertEqual(report['problems'], [])

    def test_older_python_is_rejected_even_with_tk(self):
        self.assertFalse(self.probe(version=(3, 10, 12))['ok'])

    def test_32_bit_python_is_rejected(self):
        self.assertFalse(self.probe(bits=32)['ok'])

    def test_missing_tcl_resources_are_rejected(self):
        report = self.probe(tk_error=RuntimeError('Cannot find init.tcl'))
        self.assertFalse(report['ok']);self.assertFalse(report['tkinter'])
        self.assertIn('init.tcl', report['problems'][0])

    def test_non_cpython_runtime_is_rejected(self):
        self.assertFalse(self.probe(implementation='PyPy')['ok'])

    def test_paths_with_spaces_and_cyrillic_survive_json(self):
        executable = 'C:/Пользователи/Имя с пробелами/Python311/python.exe'
        with patch.object(launch_probe.sys, 'executable', executable):
            report = self.probe()
        self.assertEqual(json.loads(json.dumps(report, ensure_ascii=True))['executable'], report['executable'])

    def test_real_probe_protocol_and_exit_status(self):
        result = subprocess.run([sys.executable, '-I', str(Path(launch_probe.__file__))], capture_output=True, text=True)
        report = json.loads(result.stdout)
        self.assertEqual(result.returncode, 0 if report['ok'] else 1)
        self.assertTrue(Path(report['executable']).is_file())


if __name__ == '__main__':
    unittest.main()
