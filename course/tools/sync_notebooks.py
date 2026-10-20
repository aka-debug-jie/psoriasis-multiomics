"""Generate paired R-kernel notebooks from the canonical Rmd lesson sources."""
from pathlib import Path
import argparse,re,textwrap,jupytext,nbformat,yaml
from nbconvert.filters.markdown_mistune import markdown2html_mistune

def render_details(source):
    # Keep each disclosure in one complete HTML block: sanitizers close isolated opening tags.
    def replace(match):
        body=markdown2html_mistune(textwrap.dedent(match.group(2)).strip())
        return '<details>'+match.group(1)+body.replace('\n','&#10;')+'</details>'
    return re.sub(r'<details>\s*(<summary>.*?</summary>)(.*?)</details>',replace,source,flags=re.S)
p=argparse.ArgumentParser();p.add_argument('--clear-outputs',action='store_true');args=p.parse_args()
root=Path(__file__).resolve().parents[1]
language=(root/'language.txt').read_text().strip() if (root/'language.txt').exists() else 'zh-CN'
for path in sorted((root/'notebooks').glob('*.Rmd')):
    notebook=jupytext.read(path)
    if notebook.cells and notebook.cells[0].cell_type=='raw' and notebook.cells[0].source.startswith('---'):
        header=yaml.safe_load(notebook.cells[0].source.strip().strip('-'))
        notebook.cells[0]=nbformat.v4.new_markdown_cell('# '+header['title'])
    notebook.metadata['kernelspec']={'name':'ir-psoriasis','display_name':'R · Multi-omics course' if language=='en' else 'R · 多组学课程','language':'R'}
    notebook.metadata['language_info']={'name':'R','file_extension':'.r','mimetype':'text/x-r','pygments_lexer':'r','codemirror_mode':'r'}
    # 两种前端各自保存练习；此脚本同步源码，不启用服务器自动配对。
    notebook.metadata.pop('jupytext',None)
    target=path.with_suffix('.ipynb')
    old=nbformat.read(target,as_version=4)if target.exists()else None
    same_code=not args.clear_outputs and old is not None and [c.source for c in notebook.cells if c.cell_type=='code']==[c.source for c in old.cells if c.cell_type=='code']
    old_codes=iter(c for c in old.cells if c.cell_type=='code')if same_code else None
    for cell in notebook.cells:
        if cell.cell_type=='markdown':cell.source=render_details(cell.source)
        if cell.cell_type=='code':
            if same_code:
                previous=next(old_codes);cell.outputs=previous.outputs;cell.execution_count=previous.execution_count
            else:cell.outputs=[];cell.execution_count=None
    nbformat.validate(notebook);nbformat.write(notebook,target)
    assert [c.source for c in jupytext.read(path).cells if c.cell_type=='code']==[c.source for c in notebook.cells if c.cell_type=='code']
print('Rmd / ipynb code pairs synchronized')
