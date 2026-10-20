# Lesson 10: enrichment, STRING networks, and drug associations

## Question and inputs

Are the seven candidate genes concentrated in particular functional terms? How are the candidate proteins connected in STRING and drug-target records? The opportunity background contains 2,045 available proteins, mapping to 2,043 unique GeneIDs. Candidates and background use the same identifiers. Each annotation family defines its own annotated background; missing annotations are not negative annotations.

## Main steps

GO annotations are propagated through `is_a` and `part_of` relations, excluding negative annotations. Reactome supplies pathway memberships. For BP, MF, CC, and Reactome, first identify the full annotated background and candidate coverage, then retain terms with 10–500 background genes. Keep terms with zero candidate overlap. The right-tail test is visible in the R code: `phyper(k-1, K, N-K, n, lower.tail=FALSE)`. N is the annotated background, K is the term's background membership, n is the annotated query, and k is the overlap. Apply BH separately within each family.

The original reference has 1,888 BP, 338 MF, 282 CC, and 342 Reactome tests, with 2, 0, 0, and 0 significant terms at FDR<0.05. The STRING plot shows associations between candidates and neighbors. Drug evidence is separated into direct targets, complex/family components, and STRING associations: 0, 5, and 53 input records respectively. Network connections do not establish direct binding or efficacy.

## Read the figures

The enrichment plot uses −log10(BH FDR), with a dashed line at FDR=0.05. Size shows candidate hits and color shows annotation family. Only terms below the cutoff are significant. Highlighted network nodes are candidates, and lines represent STRING scores; the graph is not a directional pathway. The drug plot counts distinct molecules within each class. A zero direct-target count must not be filled with network counts.

## Exercises and answers

1. **Why not use all 2,043 GeneIDs for every GO background?** Annotation coverage differs. Define N and n from each family's full annotated set before filtering term size.
2. **What changes when the upper term-size limit falls from 500 to 300?** Filter tests again and recompute their BH-adjusted values. Background coverage is fixed, and raw P values for retained terms are unchanged. Do not reuse the old FDR values.
3. **Do 53 network records mean 53 drugs directly target the candidates?** No. Direct, component, and network evidence must remain separate; network associations do not prove binding or efficacy.

Run `source("lessons/10_functional.R")`. See the [script](10_functional.R), [enrichment plot](../results/10_functional/functional_enrichment.png), [network](../results/10_functional/string_induced_network.png), and [drug-evidence plot](../results/10_functional/drug_association_classes.png).
