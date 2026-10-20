"""Build a language edition from shared statistical code and translated prose."""
from pathlib import Path
import argparse,json,re,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--language',choices=['zh-CN','en'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
source=ROOT/'course';target=a.output.resolve()
if target.exists():raise FileExistsError('Choose a new build directory; existing notebooks are never overwritten: '+str(target))
shutil.copytree(source,target,dirs_exist_ok=True,ignore=shutil.ignore_patterns('*.ipynb','*.nb.html','__pycache__','.ipynb_checkpoints','notebook_templates'))
if a.language=='en':
 overlay=ROOT/'translations/en'
 for name in ['notebooks','lessons']:
  shutil.copytree(overlay/name,target/name,dirs_exist_ok=True)
 shutil.copytree(overlay/'assets',target/'notebooks/assets',dirs_exist_ok=True)
 shutil.copy2(overlay/'COURSE.md',target/'COURSE.md')
 shutil.copytree(overlay/'data',target/'data',dirs_exist_ok=True)
translations=json.loads((ROOT/'translations/en/code_strings.json').read_text())
strings={k:v['translation'] for k,v in translations['strings'].items()}
comments={k:v['translation'] for k,v in translations['comments'].items()}
# Use the same factor labels in both editions; statistical coding is unchanged.
labels={'男':'Male','女':'Female','正常':'Normal','皮损':'Lesion'}
pattern=r'''("(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')'''
def rcode(code):
 def literal(m):
  raw=m[0];key=raw[1:-1];value=labels.get(key)
  if value is None and a.language=='en':value=strings.get(key)
  if value is None:return raw
  # Dictionary values retain R escape sequences; escape the chosen delimiter.
  return raw[0]+value.replace(raw[0],'\\'+raw[0])+raw[0]
 code=re.sub(pattern,literal,code)
 if a.language=='en':
  for old,new in sorted(comments.items(),key=lambda kv:-len(kv[0])):
   code=code.replace(old,new)
 return code
for file in target.rglob('*'):
 if not file.is_file()or file.suffix not in ('.R','.Rmd','.md','.html'):continue
 text=file.read_text(encoding='utf-8')
 if file.suffix=='.R':text=rcode(text)
 elif file.suffix in ('.Rmd','.md'):
  text=re.sub(r'(```(?:\{r[^\n]*\}|r)[^\n]*\n)(.*?)(\n```)',lambda m:m[1]+rcode(m[2])+m[3],text,flags=re.S)
 text=text.replace('课程资源','resources').replace('本课目录','chapter-folder' if a.language=='en' else '本课目录')
 file.write_text(text,encoding='utf-8')
(target/'language.txt').write_text(a.language+'\n')
import csv,hashlib
with(target/'data/manifest.tsv').open('w',newline='',encoding='utf-8')as handle:
 writer=csv.DictWriter(handle,['path','bytes','sha256'],delimiter='\t');writer.writeheader()
 for file in sorted((target/'data').rglob('*')):
  if file.is_file()and file.name!='manifest.tsv':writer.writerow({'path':file.relative_to(target).as_posix(),'bytes':file.stat().st_size,'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
subprocess.run([sys.executable,str(target/'tools/sync_notebooks.py'),'--clear-outputs'],check=True)
print('Prepared',a.language,target)
