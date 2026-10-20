# External tools

Student releases contain native SMR and PLINK executables. Their actual paths,
versions and hashes are written into `manifest.tsv` when the package is assembled;
no tool is downloaded during student startup.

SMR is built from fixed commit `eedb1a62ba4c66ab0b09acd0205876955a803f18`, with the
course's exact-|Z|-tie ordering patch and Windows I/O changes where applicable.
See [the method note](../notes/smr-portability.md). A stock, unpatched SMR may choose
different representative SNPs and is not the course's reference build.

The source-repository manifest describes the required tools but does not include
executables. Maintainer builds use `scripts/build_native.py` from the repository
root. For a separately configured source checkout, `TEACHING_SMR` and
`TEACHING_PLINK` can point to matching local binaries. The upstream SMR GPL version
2 file and PLINK GPL version 3 file are retained here; further dependency notices
and corresponding source accompany the releases.
