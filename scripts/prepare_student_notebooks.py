"""Prepare empty student workbooks after exporting the executed HTML pages."""
from pathlib import Path
import argparse
import nbformat


def prepare(course):
    templates = course / 'notebook_templates'
    templates.mkdir(exist_ok=True)
    count = 0
    for source in sorted((course / 'notebooks').glob('[0-9][0-9]_*.Rmd')):
        path = source.with_suffix('.ipynb')
        if not source.with_suffix('.nb.html').is_file():
            raise RuntimeError(f'Export the executed reading page first: {source.name}')
        nb = nbformat.read(path, as_version=4)
        nb.metadata.pop('widgets', None)
        for cell in nb.cells:
            if cell.cell_type == 'code':
                cell.outputs = []
                cell.execution_count = None
                cell.metadata.pop('execution', None)
        nbformat.validate(nb)
        nbformat.write(nb, path)
        nbformat.write(nb, templates / path.name)
        count += 1
    return count


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('course', type=Path)
    print('Prepared empty workbooks:', prepare(parser.parse_args().course))
