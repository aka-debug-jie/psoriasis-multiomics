"""Add diagnostic output only; never used for student release builds."""
from pathlib import Path
import sys

path = Path(sys.argv[1]) / 'src/SMR_data.cpp'
source = path.read_text()
start = source.index('double heidi_test_new(bInfo*')
end = source.index('int get_num_probes(', start)
section = source[start:end]
anchor = '  nsnp = sn_ids.size();\n  double pdev = -9;'
assert section.count(anchor) == 1
trace = r'''
  for (int j = 0; j < sn_ids.size(); ++j) {
    const int k = sn_ids[j];
    fprintf(stderr, "TRACE_INPUT\t%d\t%d\t%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
            (int)smrwk_heidi.cur_prbidx, j, smrwk_heidi.rs[k].c_str(),
            smrwk_heidi.zxz[k], _bxz[j], _sexz[j], _byz[j], _seyz[j]);
    for (int i = 0; i < sn_ids.size(); ++i)
      fprintf(stderr, "TRACE_LD\t%d\t%d\t%d\t%.17g\n",
              (int)smrwk_heidi.cur_prbidx, j, i, C(j, i));
  }
'''
source = source[:start] + section.replace(anchor, trace + anchor) + source[end:]
start = source.index('void smr_heidi_func(std::vector<SMRRLT>& smrrlts, char* outFileName, bInfo*')
end = source.index('void smr_heidi_func(', start + 1)
section = source[start:end]
anchor = '      long nsnp = -9;'
assert section.count(anchor) == 1
section = section.replace(anchor, '      fprintf(stderr, "TRACE_PROBE\\t%d\\t%s\\n", (int)smrwk.cur_prbidx, probename.c_str());\n' + anchor)
source = source[:start] + section + source[end:]
path.write_text(source)
print('Added HEIDI input and LD diagnostics without changing selection or statistics.')
