# Ordering tied instruments in SMR

The course uses a fixed SMR 1.4.3 source revision with one explicit rule for HEIDI instrument selection: exactly equal `|Z|` values are ordered by their input indices within the current selection step. LD limits, instrument thresholds, and statistical formulas retain their specified values.

## Why the order matters

In the ERAP2 example, rs12189125 and rs1974871 have identical expression effects and standard errors, with `Z=-28.7809042484` and reference `r²=0.99237`. Another exact tie, rs1230361 and rs1235617, has identical reference genotypes. Their disease-GWAS standard errors differ.

The original comparator uses only `|Z|`, leaving tied order unspecified. Standard libraries can return different orders. Subsequent LD pruning uses indices to choose between nodes of equal degree, so the retained representative variant can differ.

In this comparison, the original builds entered pruning with the same 352 variants, 7,841 correlation edges, and node degrees. The tied index order exactly explained the differing removal lists. The course comparator is:

```cpp
const double za = fabs(smrwk->zxz[a]), zb = fabs(smrwk->zxz[b]);
return za == zb ? a < b : za > zb;
```

## Comparison on the default inputs

| Tool and environment | ERAP2 HEIDI P | LIPA HEIDI P |
|---|---:|---:|
| Original, Linux / Windows | 0.2480520 | 0.2011178 |
| Original, Intel / Apple silicon Mac | 0.2465839 | 0.2011178 |
| Course build, all four environments | 0.2423640 | 0.2011178 |

Each environment was run twice with one thread and twice with two threads. All reported fields from the course builds matched. Windows uses CRLF line endings; after normalizing line endings, the complete file hashes matched as well. HEIDI retained 20 variants for ERAP2 and 9 for LIPA. SMR effects, standard errors, and P values retained their original values.

The classroom reference uses the course build. The original Linux result is kept separately in `data/smr/reference_original_linux.tsv`. The research snapshot retains its original tools and analysis records. This rule makes the selection order explicit; expression association and heterogeneity still require separate interpretation.

The base source revision is `eedb1a62ba4c66ab0b09acd0205876955a803f18`. The release supplies the patch, actual source, and build settings. In the project repository, `scripts/patch_smr_order.py` specifies the change, and `smr-portability.yml` compares the original, diagnostic-only, and explicitly ordered builds.
