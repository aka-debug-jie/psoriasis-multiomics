# Public course inputs

These are the actual inputs for course v4. `manifest.tsv` records file identities,
and `sources.tsv` records sources. Provider-specific licences and attribution
requirements are separate from the MIT licence of the original course code.

| Folder | Data and starting point |
|---|---|
| bulk | NCBI-reprocessed integer counts for GSE54456, GRCh38.p13; 171 samples (90 lesions, 81 normal). These are neither the submitter's RPKM nor the paper's historical count matrix. |
| genetics | Public exposure records for a fixed family of 2,709 eligible proteins. MVP is the discovery outcome and Stuart supplies external validation. MVP SE is reconstructed from the reported CI; the source P is retained separately. Raw FinnGen subsets are not distributed in this edition. |
| smr | eQTL Catalogue reprocessed GTEx blood QTD000356: ERAP2 and LIPA teaching examples, 670 participants and variant-level an/2 of 636–670. ESD, fixed BESD, MVP regional data and 1000G EUR LD are provided. Genes were chosen using exposure eQTL strength, not as new discovered targets. |
| coloc | Fenland Olink IL7R protein and GTEx ERAP2 mRNA are separately compared with MVP. Inputs retain available, shared, harmonizable regional SNPs. A requested 1 Mb window can exceed source cis coverage; see regions.tsv. |
| single_cell | GSE173706, 22 samples, at most 500 cells per sample sampled before QC with seed=20261017; 10,708 cells. Pseudobulk inputs aggregate all 46,225 QC-passing research cells by sample and cell type and do not change with classroom reclustering. |
| pgs | 503 EUR reference individuals from 1000 Genomes Phase3: two IL7R score SNPs and a separate set of 5,000 PCA SNPs. Weights come from one UKB-PPP protein study. These reference individuals have no disease follow-up for this course. |
| cox | 228 NCCTG lung-cancer observations from survival::lung. Time is in days; status 1 is censoring and 2 is death. Missing values remain missing. These are not the individuals in the PGS chapter. |
| functional | Seven research candidates and an opportunity background of 2,043 GeneIDs; GO and Reactome release 92 memberships. Reactome replaces the undistributed research KEGG snapshot. STRING associations and direct drug-target records from ChEMBL remain distinct. |
| clinical | Fixed instrument-position statistics for 22 MVP European-ancestry PheCode outcomes. Accessions, exact phenotype names, case/control counts and differences are in endpoints.tsv. This is a separate teaching panel, not a phenotype-equivalent replacement of the research FinnGen cohort. |
| evidence | Author-reported score flags for 29 targets and a separate snapshot of the seven-candidate research results. They teach evidence integration and are not automatically updated by the smaller public examples. Missing UKB evidence remains NA. |

SMR coordinates are uniquely mapped from GRCh38 to GRCh37, complementing alleles
on negative-strand mappings and excluding ambiguous mappings or duplicate IDs.
Colocalization further excludes unresolved palindromic or incompatible alleles.
Source metadata, mapping rules and input identities are recorded in JSON files.

MVP outcomes can share participants, controls, and GWAS information with discovery.
Several supported cells in a matrix are not several independent replications.
Research result snapshots and public teaching inputs are distinguished explicitly.
Raw FinnGen/eQTLGen subsets and the KEGG membership snapshot are excluded from
the public package because redistribution was not established in this review.
