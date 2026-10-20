"""Read upstream sources from recorded recipes, including rattler staging caches."""
from pathlib import Path
import yaml


def recipe_sources(directory):
    for name in ('rendered_recipe.yaml', 'meta.yaml', 'recipe.yaml'):
        path = Path(directory) / name
        if not path.is_file():
            continue
        try:
            document = yaml.safe_load(path.read_text())
        except yaml.YAMLError:
            continue
        if not isinstance(document, dict):
            continue
        recipe = document.get('recipe', document) if name == 'rendered_recipe.yaml' else document
        entries = []
        for section in [recipe, *recipe.get('staging_caches', [])]:
            source = section.get('source', [])
            entries.extend([source] if isinstance(source, dict) else source)
        if entries:
            return entries, path
    return [], None
