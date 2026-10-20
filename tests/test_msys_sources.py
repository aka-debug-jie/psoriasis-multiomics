"""A locally built DLL must not be attributed to an unrelated system package."""
from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch

from scripts.collect_msys_sources import collect_msys_sources


class DllOriginTests(unittest.TestCase):
    def test_smr_dependency_requires_matching_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            dlls, deps, notices, sources = [root / name for name in
                                           ('bin', '_deps', 'notices', 'sources')]
            for folder in (dlls, deps / 'zlib-build', notices, sources):
                folder.mkdir(parents=True)
            (dlls / 'libzlib.dll').write_bytes(b'compiled zlib')
            built = deps / 'zlib-build/libzlib.dll'
            built.write_bytes(b'compiled zlib')
            with patch('scripts.collect_msys_sources.subprocess.check_output') as query:
                collect_msys_sources(root / 'msys', dlls, notices, sources, deps)
                query.assert_not_called()
            record = json.loads((notices / 'msys2-components.json').read_text())['libraries'][0]
            self.assertEqual(record['owner'], 'SMR build dependency')
            self.assertEqual(record['source_archive'], 'SMR-1.4.3-complete-source.tar.gz')
            built.write_bytes(b'a different build')
            with self.assertRaisesRegex(AssertionError, 'Unresolved DLL origin'):
                collect_msys_sources(root / 'msys', dlls, notices, sources, deps)


if __name__ == '__main__':
    unittest.main()
