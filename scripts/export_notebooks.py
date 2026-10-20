"""Create offline reading pages from genuinely executed R notebooks."""
from pathlib import Path
import argparse,nbformat
from nbconvert import HTMLExporter
from bs4 import BeautifulSoup
p=argparse.ArgumentParser();p.add_argument('course',type=Path);p.add_argument('--language',choices=['zh-CN','en'],default='zh-CN');p.add_argument('--lesson',choices=[f'{i:02d}'for i in range(13)]);a=p.parse_args()
exporter=HTMLExporter(template_name='classic',embed_images=True)
exporter.exclude_input_prompt=True;exporter.exclude_output_prompt=True
for file in sorted((a.course/'notebooks').glob((a.lesson+'_'if a.lesson else'')+'*.ipynb')):
 nb=nbformat.read(file,as_version=4)
 assert all(c.execution_count is not None for c in nb.cells if c.cell_type=='code'),file
 assert not any(o.output_type=='error'for c in nb.cells for o in c.get('outputs',[])),file
 body,_=exporter.from_notebook_node(nb,resources={'metadata':{'path':str(file.parent)}})
 soup=BeautifulSoup(body,'html.parser')
 # R text/table/PNG outputs are static. They need no remote widget, math, or JS loader.
 for node in soup.find_all('script'):node.decompose()
 for node in soup.find_all('link'):
  if node.get('href','').startswith(('http:','https:','//')):node.decompose()
 for node in soup.find_all(['img','iframe','source']):
  if node.get('src','').startswith(('http:','https:','//')):raise RuntimeError('Remote display asset: '+node['src'])
 for node in soup.find_all('a'):
  if node.get('href')=='00_start.ipynb' :node['href']='00_start.nb.html'
 soup.html['lang']=a.language
 title=nb.cells[0].source.lstrip('# ').splitlines()[0]
 if soup.title:soup.title.string=title
 banner=soup.new_tag('p');banner['style']='padding:12px 16px;background:#eaf4ef;border-radius:8px;color:#185e4b'
 banner.string='Executed course example · These outputs use the default parameters. Open a blank notebook in Jupyter to run the code yourself.' if a.language=='en'else'已运行的课程展示 · 以下输出来自默认参数示例。亲手练习请打开 Jupyter 中的空白 Notebook。'
 container=soup.find(id='notebook-container')
 if container is not None:container.insert(0,banner)
 style=soup.new_tag('style');style.string='''body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Microsoft YaHei",sans-serif;color:#23342f;background:#f5f8f6}#notebook-container{max-width:1120px;box-shadow:none;border:1px solid #dce7e1;border-radius:12px;padding:30px}.rendered_html{font-size:16px;line-height:1.85}.rendered_html h1,.rendered_html h2,.rendered_html h3{color:#185e4b}.input_area{border-radius:6px!important;background:#f4f7f6!important}.output_area{margin-bottom:12px}pre{overflow:auto}.text_cell_render{padding:12px 6px}a{color:#17684f}details{margin:12px 0;padding:10px 14px;background:#f4f8f5;border-left:3px solid #7ca992;border-radius:6px}summary{cursor:pointer;color:#185e4b;font-weight:600}details[open]>summary{margin-bottom:10px}.rendered_html table{display:block;max-width:100%;overflow-x:auto}@media(max-width:600px){#notebook-container{padding:12px!important}.rendered_html{font-size:15px}.input_area{overflow-x:auto}}''';soup.head.append(style)
 file.with_suffix('.nb.html').write_text(str(soup),encoding='utf-8')
print('Offline displays exported')
