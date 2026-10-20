"""Build-time output stripping must preserve the independently exported reading page."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import shutil
import subprocess
import sys
import nbformat
from nbconvert import HTMLExporter
from bs4 import BeautifulSoup

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/prepare_student_notebooks.py'
spec = importlib.util.spec_from_file_location('prepare_student_notebooks', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class StudentNotebookTests(unittest.TestCase):
    def test_disclosure_body_survives_notebook_and_html_conversion(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'notebooks').mkdir();(root/'tools').mkdir()
            source=Path(__file__).resolve().parents[1]/'course/tools/sync_notebooks.py'
            shutil.copy2(source,root/'tools/sync_notebooks.py')
            (root/'notebooks/00_start.Rmd').write_text(
                '# Read code\n\n<details><summary>Explain</summary>\n\n'
                '`match()` finds a row.\n\n| Name | Value |\n|---|---|\n| A | 1 |\n\n'
                '</details>\n\n```{r}\n1 + 1\n```\n')
            subprocess.run([sys.executable,str(root/'tools/sync_notebooks.py')],check=True,capture_output=True)
            nb=nbformat.read(root/'notebooks/00_start.ipynb',as_version=4)
            self.assertEqual([c.source for c in nb.cells if c.cell_type=='code'],['1 + 1'])
            html,_=HTMLExporter(template_name='classic').from_notebook_node(nb)
            soup=BeautifulSoup(html,'html.parser');details=soup.find('details')
            self.assertIn('finds a row',details.get_text())
            self.assertIsNotNone(details.find('table'))
            self.assertIsNone(details.find(class_='input_area'))
            self.assertFalse(details.has_attr('open'))

    def test_reading_output_and_student_source_are_separate(self):
        with tempfile.TemporaryDirectory() as directory:
            course = Path(directory)
            notebooks = course / 'notebooks'
            notebooks.mkdir()
            source = notebooks / '00_start.Rmd'
            source.write_text('```{r}\n1 + 1\n```\n')
            html = source.with_suffix('.nb.html')
            html.write_text('<p>Real execution result: 2</p>')
            nb = nbformat.v4.new_notebook(cells=[
                nbformat.v4.new_markdown_cell('Run this cell.'),
                nbformat.v4.new_code_cell('1 + 1', execution_count=1,
                    outputs=[nbformat.v4.new_output('stream', name='stdout', text='2')],
                    metadata={'execution': {'iopub.status.busy': 'timestamp'}})
            ])
            nb.metadata['widgets'] = {'old': 'state'}
            nbformat.write(nb, source.with_suffix('.ipynb'))
            self.assertEqual(module.prepare(course), 1)
            actual = nbformat.read(source.with_suffix('.ipynb'), as_version=4)
            self.assertEqual(actual.cells[1].source, '1 + 1')
            self.assertEqual(actual.cells[1].outputs, [])
            self.assertIsNone(actual.cells[1].execution_count)
            self.assertNotIn('execution', actual.cells[1].metadata)
            self.assertNotIn('widgets', actual.metadata)
            self.assertEqual(html.read_text(), '<p>Real execution result: 2</p>')
            self.assertEqual(source.with_suffix('.ipynb').read_bytes(),
                (course / 'notebook_templates/00_start.ipynb').read_bytes())

    def test_export_is_required_before_stripping(self):
        with tempfile.TemporaryDirectory() as directory:
            course = Path(directory)
            (course / 'notebooks').mkdir()
            (course / 'notebooks/00_start.Rmd').write_text('')
            with self.assertRaisesRegex(RuntimeError, 'Export the executed'):
                module.prepare(course)


if __name__ == '__main__':
    unittest.main()
