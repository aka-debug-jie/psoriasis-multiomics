"""Run the unchanged two-probe example repeatedly, separately from course packaging."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--smr', type=Path, required=True)
parser.add_argument('--platform', required=True)
parser.add_argument('--label', choices=['baseline', 'trace'], required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
data = root / 'course/data/smr'
out = root / '.build/smr-diagnostics' / args.platform / args.label
out.mkdir(parents=True, exist_ok=True)
smr = args.smr.resolve()
inputs = ['EUR_two_regions.bed', 'EUR_two_regions.bim', 'EUR_two_regions.fam',
          'public_eqtl.besd', 'public_eqtl.epi', 'public_eqtl.esi', 'mvp_regional.ma', 'probes.txt']
record = {'platform': platform.platform(), 'label': args.label,
          'binary_sha256': hashlib.sha256(smr.read_bytes()).hexdigest(),
          'inputs': {name: hashlib.sha256((data / name).read_bytes()).hexdigest() for name in inputs},
          'status': 'running', 'runs': []}
(out / 'diagnostic.json').write_text(json.dumps(record, indent=2) + '\n')
def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

for threads in (1, 2):
    for repeat in (1, 2):
        prefix = out / f'threads-{threads}-repeat-{repeat}'
        command = [str(smr), '--bfile', str(data / 'EUR_two_regions'),
                   '--gwas-summary', str(data / 'mvp_regional.ma'),
                   '--beqtl-summary', str(data / 'public_eqtl'),
                   '--peqtl-smr', '5e-8', '--heidi-skip-thresh', '1',
                   '--peqtl-heidi', '0.0015654', '--ld-upper-limit', '0.9',
                   '--ld-lower-limit', '0.05', '--cis-wind', '2000',
                   '--diff-freq', '0.2', '--diff-freq-prop', '0.05',
                   '--thread-num', str(threads), '--extract-probe', str(data / 'probes.txt'),
                   '--out', str(prefix)]
        with prefix.with_suffix('.log').open('w') as handle:
            subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
        rows = read_table(prefix.with_suffix('.smr'))
        record['runs'].append({'threads': threads, 'repeat': repeat, 'results': rows})
        (out / 'diagnostic.json').write_text(json.dumps(record, indent=2) + '\n')
        if repeat == 2:
            assert rows == read_table(out / f'threads-{threads}-repeat-1.smr'), 'Repeat differs'
        if args.label == 'trace':
            assert rows == read_table(out.parent / 'baseline' / prefix.with_suffix('.smr').name), 'Tracing changed output'
        print(args.label, threads, repeat, [(x['Gene'], x['p_HEIDI'], x['nsnp_HEIDI']) for x in rows], flush=True)
record['status'] = 'passed'
(out / 'diagnostic.json').write_text(json.dumps(record, indent=2) + '\n')
