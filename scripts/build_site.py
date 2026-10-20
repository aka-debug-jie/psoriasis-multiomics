"""Build a static bilingual course site from validated, executed reading pages."""
from pathlib import Path
import argparse,json,shutil,html,re
from bs4 import BeautifulSoup
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--courses',type=Path,required=True);p.add_argument('--output',type=Path,default=ROOT/'site-output');a=p.parse_args()
release=json.loads((ROOT/'config/release.json').read_text());a.output.mkdir(parents=True,exist_ok=True)
(a.output/'.nojekyll').touch()
(a.output/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Psoriasis multi-omics</title><meta http-equiv="refresh" content="0;url=en/"><p><a href="en/">English</a> · <a href="zh/">简体中文</a></p></html>')
repo='https://github.com/'+release['repository']
for language,slug in [('en','en'),('zh-CN','zh')]:
 course=a.courses/language;target=a.output/slug;target.mkdir(exist_ok=True)
 destination=target/'course'
 shutil.copytree(course,destination,dirs_exist_ok=True,ignore=shutil.ignore_patterns('*.log','*.rds','*.RData','__pycache__','.ipynb_checkpoints','notebook_templates','tools'))
 for file in destination.rglob('*.nb.html'):
  s=file.read_text();s=re.sub(r'href="([0-9]{2}_[^"/]+)\.ipynb"',r'href="\1.nb.html"',s);file.write_text(s)
 shutil.copy2(ROOT/'site/style.css',target/'style.css')
 en=language=='en'
 title='Learn multi-omics, one R cell at a time'if en else'逐格运行 R，学会多组学分析'
 intro=('<p>A practical course built around psoriasis research. Read the worked examples, run the code, and learn to explain what each result supports.</p>'if en else'<p>以银屑病研究为线索，从读懂输入、逐格运行代码，到解释统计结果。先看已运行的示例，再亲手完成分析与参数练习。</p>')
 nav=''.join(f'<a class="entry-card" href="{url}"><strong>{label}</strong><span>{note}</span></a>'for url,label,note in [
  ('course/notebooks/00_start.nb.html','Start reading'if en else'从第 0 章开始','Worked examples, with code and outputs'if en else'讲解、代码和实际运行结果'),
  ('#downloads-title','Run the course'if en else'下载动手练习','Bundled R and Jupyter for your computer'if en else'自带 R 和 Jupyter 的课程包'),
  (repo+'/tree/main/research','Research sources'if en else'查看科研方法','Code, provenance, results, and differences'if en else'代码、来源、结果和与原文的差异')])
 cards=[]
 for file in sorted((destination/'notebooks').glob('[0-9][0-9]_*.nb.html')):
  soup=BeautifulSoup(file.read_text(),'html.parser');heading=soup.find('h1')
  label=heading.get_text(' ',strip=True)if heading else file.stem
  label=re.sub(r'^[0-9]{2}[\s｜|:：]*','',label)
  cards.append(f'<a class="course-card" href="course/notebooks/{file.name}"><span class="chapter-number">{file.name[:2]}</span><h3>{html.escape(label)}</h3></a>')
 downloads=[]
 for platform,label in [('windows-x64','Windows x64'),('macos-intel','Mac · Intel'),('macos-arm64','Mac · Apple Silicon')]:
  assets=[x for x in release['assets']if x['platform']==platform and x['language']==language and x.get('status')=='ready']
  if assets:
   asset=assets[0];url=repo+'/releases/download/v'+release['version']+'/'+asset['name']
   button=f'<a class="download-button" href="{url}">{"Download ZIP"if en else"下载 ZIP"}</a>'
   note=f"v{release['version']} · {asset['bytes']/1e6:.1f} MB · "+('Native tests passed'if en else'原生运行检查通过')
   checksum=f'<details><summary>SHA256</summary><code>{asset["sha256"]}</code></details>'
  else:
   button=f'<span class="download-button is-disabled" aria-disabled="true">{"Preparing release"if en else"正在准备发布"}</span>'
   note='Download will open after native validation.'if en else'完成该系统的原生运行验收后开放下载。';checksum=''
  downloads.append(f'<article class="download-card"><div><h3>{label}</h3><p>{note}</p>{checksum}</div>{button}</article>')
 extra=('<p>Extract the complete ZIP, then open <code>Run course.bat</code> on Windows or <code>Run course.command</code> on Mac. The package includes R 4.4.2 and Jupyter. Choose the correct Mac chip; check About This Mac if unsure. Use <code>View course.html</code> to read without starting R.</p>'if en else'<p>完整解压 ZIP 后，Windows 打开 <code>Run course.bat</code>，Mac 打开 <code>Run course.command</code>。课程已含 R 4.4.2 和 Jupyter。Mac 请按“关于本机”中的芯片选择版本。只阅读时打开 <code>View course.html</code>。</p>')
 research=('<p>The research code independently reconstructs the analyses in <a href="https://doi.org/10.1093/bioinformatics/btag347">btag347</a>. Its seven-candidate result uses documented substitute inputs and method adjustments; it is not an exact reproduction of the paper’s 29 candidates. UK Biobank individual-level PGS and survival analyses were not reproduced.</p><p>The public course uses separately documented teaching inputs where needed. The 1000 Genomes scoring exercise and the survival::lung example are independent datasets. Source versions, licences, and changes are recorded with the course.</p>'if en else'<p>科研代码独立重建了 <a href="https://doi.org/10.1093/bioinformatics/btag347">btag347</a> 的分析流程。七个候选来自已记录的替代输入和方法修正，与论文的 29 个候选分别对照；UK Biobank 个体层面的遗传评分和生存分析未复现。</p><p>公开课程在必要环节采用另行说明的教学输入。1000 Genomes 遗传评分和 survival::lung 生存分析使用独立数据，不能拼成同一队列。数据版本、许可和变化均随课程记录。</p>')
 footer=f'<p><a href="{repo}">GitHub</a> · <a href="{repo}/blob/main/LICENSE">MIT code</a> · <a href="{repo}/blob/main/LICENSE-CONTENT.md">CC BY 4.0 content</a> · <a href="{repo}/blob/main/DATA_LICENSES.tsv">Data sources & licences</a></p>'
 source=(ROOT/'site/home.html').read_text()
 for key,value in {'LANG':language,'TITLE':title,'INTRO':intro,'NAV':nav,'COURSE_CARDS':''.join(cards),'DOWNLOADS':''.join(downloads)+extra,'RESEARCH':research,'FOOTER':footer}.items():source=source.replace('{{'+key+'}}',value)
 soup=BeautifulSoup(source,'html.parser')
 for node in list(soup.find_all('span',lang=True)):
  if node.get('lang')!=language:node.decompose()
 for node in soup.find_all('span',attrs={'aria-hidden':'true'}):
  if node.get_text(strip=True)=='/':node.decompose()
 for link in soup.select('.language-nav a'):link['href']='../en/'if link.get('lang')=='en'else'../zh/'
 (target/'index.html').write_text(str(soup))
size=sum(f.stat().st_size for f in a.output.rglob('*')if f.is_file())
assert size<950*1024**2, 'Static site is too large for Pages'
print('Static bilingual site:',size,'bytes')
