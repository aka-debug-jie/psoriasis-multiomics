# Preparing the public teaching inputs

Student packages contain the reviewed data snapshot. These commands are for maintainers rebuilding inputs; students do not run them.

The release snapshot is identified by `config/data_asset.json`, and its members by `course/data/manifest.tsv`. Source identities and reuse terms are listed in `course/data/sources.tsv` and `DATA_LICENSES.tsv`. The research input pipeline remains separate under `research/`.

## GTEx whole-blood examples

Use Python with requests, pysam, pandas and pyliftover. The pinned Ensembl metadata fixes the gene coordinates used in the released example:

```sh
python scripts/fetch_public_eqtl.py --output .build/eqtl-inputs --metadata-dir config/eqtl_metadata
python scripts/prepare_public_eqtl.py --course course --records .build/eqtl-inputs --chain /path/to/hg38ToHg19.over.chain.gz --mvp /path/to/mvp.full.ma --ld /path/to/1000G_EUR_prefix --plink /path/to/plink --smr /path/to/smr
```

The downloader reads indexed regions from QTS000015/QTD000356 rather than the complete study file. Its receipt records URLs and hashes. In release preparation, it recovered exactly 7,999 ERAP2 and 8,027 LIPA records, with byte-identical TSVs to the input snapshot used for mapping. ERAP2 and LIPA were selected using exposure strength; they are not additional discovery candidates.

The preparation step uses a unique GRCh38-to-GRCh37 mapping, records exclusions, complements alleles on the negative strand and prepares ESD/BESD, MVP regional and LD inputs. Colocalization excludes ambiguous palindromic or incompatible alleles. The source cis coverage does not fully span the requested ERAP2 1 Mb interval; `regions.tsv` retains the actual boundary.

## Other replacement inputs

`python scripts/prepare_reactome.py --help` describes the fixed Reactome release 92 preparation. `python scripts/prepare_public_clinical.py --help` describes acquisition of the 22 MVP European-ancestry outcomes in `config/clinical_public_panel.tsv`. These retain the study phenotypes and sample counts, rather than treating the panel as a one-to-one FinnGen replacement.

After changing inputs, rerun the affected analyses, rebuild the manifest and regenerate both language editions. The fixed evidence-integration snapshot in Chapter 12 is not automatically overwritten by teaching examples.
