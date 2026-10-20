"""Resolve exactly tied |Z| values by input index; thresholds are unchanged."""
from pathlib import Path
import sys

path = Path(sys.argv[1]) / 'src/SMR_data.cpp'
source = path.read_text()
old = 'auto cmp = [&](int a, int b) { return fabs(smrwk->zxz[a]) > fabs(smrwk->zxz[b]); };'
new = '''auto cmp = [&](int a, int b) {
      const double za = fabs(smrwk->zxz[a]), zb = fabs(smrwk->zxz[b]);
      return za == zb ? a < b : za > zb;
    };'''
assert source.count(old) == 1, 'Expected the pinned SMR comparator exactly once'
notice = '// Modified 2026-10-21 by the psoriasis-multiomics course: deterministic input-index ordering for exact |Z| ties.\n'
path.write_text(notice + source.replace(old, new))
print('Exactly tied |Z| values now follow their original input indices.')
