from pathlib import Path
import tempfile
import unittest

from scripts.conda_recipe_sources import recipe_sources


class RecipeSourceTests(unittest.TestCase):
    def test_rendered_split_package_includes_staging_source(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'recipe.yaml').write_text('source:\n  url: https://example.org/${{ version }}.tar.gz\n')
            (root / 'rendered_recipe.yaml').write_text(
                'recipe:\n  staging_caches:\n  - name: glib-build\n'
                '    source:\n    - url: https://example.org/glib-2.90.1.tar.xz\n'
                '      sha256: abc\nfinalized_dependencies:\n'
                '  url: https://example.org/not-a-source.conda\n')
            entries, path = recipe_sources(root)
            self.assertEqual(path.name, 'rendered_recipe.yaml')
            self.assertEqual(entries, [{'url': 'https://example.org/glib-2.90.1.tar.xz', 'sha256': 'abc'}])

    def test_classic_source_dictionary(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'meta.yaml'
            path.write_text('source:\n  url: https://example.org/source.tar.gz\n')
            self.assertEqual(recipe_sources(path.parent)[0], [{'url': 'https://example.org/source.tar.gz'}])


if __name__ == '__main__':
    unittest.main()
