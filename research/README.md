# Research reconstruction

This directory contains the analysis code and an aggregate result snapshot for
an independent reconstruction of *Integrated multi-omics strategies for
identifying novel therapies in psoriasis* (Bioinformatics, btag347).

The recorded primary analysis tested a fixed family of 2,709 proteins and produced
seven candidates: IL12B, SGSH, TNFAIP3, DDX58, IL7R, PRSS53 and UBLCP1. It accounts
for repeated use of the same disease GWAS record in the protein-source combination.
Sensitivity analyses, expression results and evidence coverage are reported
separately. These are reconstruction results under documented source and method
choices, not the paper's original 29-candidate result.

Start with the [result overview](final/README.md), [research report](final/REPORT.md),
[figure atlas](final/FIGURE_ATLAS.pdf), [table dictionary](final/DATA_DICTIONARY.md)
and [source list](final/SOURCES.tsv). The report is currently in Chinese; table
fields, code, and this overview provide an English entry point. The teaching
course has a complete separate English edition.

## Data and scope

The research snapshot preserves the original project choices, including FinnGen
and eQTLGen analyses. Raw copies of those datasets are not redistributed here.
The public teaching edition uses independently documented examples where
redistribution was not established. It does not overwrite this research snapshot.

Important remaining differences are the substituted Fenland Olink input, unavailable
deCODE full regional statistics, 171/174 bulk samples, reconstructed single-cell
parameters and broad cell labels, current database snapshots, and the unresolved
MVP source-P versus CI-derived-SE convention. UK Biobank individual-level PGS/Cox
analyses were not reproduced. The course's 1000 Genomes and lung examples do not
fill those research evidence slots.

## Running research code

The research scripts target the recorded Linux environment and require the
originally documented source inputs. Environment folders, large raw downloads,
restricted materials and working logs are excluded from Git. Student packages
have a separate portable environment and do not require these steps.

From the repository root:

```bash
bash research/scripts/bootstrap_env.sh
bash research/scripts/bootstrap_data_env.sh
bash research/scripts/bootstrap_sc_env.sh
research/envs/data/bin/python research/scripts/prepare_mr_methods.py
research/envs/r442_sc/bin/Rscript --vanilla research/execute.R --list
```

Use the source and download records in `config/` and the acquisition scripts to
prepare the selected module inputs under `research/data/`, respecting each
provider's terms. Then, for example:

```bash
research/envs/r442_sc/bin/Rscript --vanilla research/execute.R --module mr_substitute
research/envs/r442_sc/bin/Rscript --vanilla research/execute.R --module final_research --resume
```

The integration step consumes upstream module outputs; it does not silently rerun
all underlying analyses. `--resume` checks recorded input/output/script identities.
`final/PUBLIC_MANIFEST.json` identifies the published aggregate files. The original
full-input validation and published-snapshot checks have different scopes.
