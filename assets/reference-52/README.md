# 52-sample reference figures

[View the gallery in the project README](../../README.md#52-sample-reference-gallery).

These 18 figures come from the integrated 2026-11-03 analysis of
52 *Geobacillus*/*Parageobacillus* genomes (40 and 12, respectively).
The source DNMBcluster revision was
[`7d19175`](https://github.com/JAEYOONSUNG/DNMBcluster/commit/7d191757b0829553a99930a9612e3e21cf81c3bd).

## Run and dataset

| Item | Value |
|---|---|
| Genomes | 52 |
| Annotated proteins | 177,974 |
| Pan-genome clusters | 14,337 |
| Core clusters | 1,034 |
| Engine / sequence level | MMseqs2 / protein |
| Identity / coverage thresholds | 0.4 / 0.8 |
| Alignment metrics | Bidirectional alignment enabled |

Clustering settings from the recorded run, with portable input/output paths:

```bash
dnmbcluster run ./genbanks_52 \
  --output ./results_52 \
  --tool mmseqs2 --level protein \
  --identity 0.4 --coverage 0.8 \
  --threads 10 --max-ram 24G --alignment
```

The figures include later display refinements to the flower plot and MDS
(strain names, species colors, and layout), and supplementary similarity
matrices in the Circos plot. Re-running the command above alone does not
reproduce those refinements or the supplementary alignment analysis.

## Supplementary ANIb and POCP

Both symmetric matrices cover all 1,326 unordered genome pairs.

- **ANIb:** BLASTn 2.17.0+ with pyANI 0.2.13.1's 1,020-nt fragmentation
  and filtering; the displayed similarity is the arithmetic mean of the
  two directional estimates.
- **POCP:** direct BLASTP 2.17.0+ against complete proteomes. A protein is
  conserved when a reported HSP has E-value < 1e-5, identity > 40%, and
  alignment length / query protein length > 0.50, with SEG disabled.
  `POCP(A,B) = 100 * (C_A_to_B + C_B_to_A) / (T_A + T_B)`.

These are supplementary alignment results, distinct from DNMBcluster's
native `--ani` stage (skani ANI and shared-cluster POCP). Circos genome
ordering uses average-linkage clustering on `100 - ANIb`; the ordering
is a similarity dendrogram, not an inferred phylogenetic tree.

## Files and provenance

Original PDFs are copied byte-for-byte from the final integrated
52-genome figure set. Existing PNGs for the flower plot, MDS, Circos,
and the two similarity heatmaps are also copied unchanged. Other PNGs
are rendered from the PDFs with Poppler `pdftoppm`, using 2,400 pixels
on the longest side. The UpSet PDF has a blank first page; its preview
uses page 2. Original plot labels and layout are preserved.

[manifest.json](manifest.json) records source figure names, preview pages,
conversion methods, and SHA-256 hashes for every included PDF and PNG.
Run details were checked against the source package's `PROTOCOL.md`,
`canonical_input_manifest.tsv`, `figure_validation.json`, and
`initial_stage_audit/missing_stage_report.json`.

| Figure | Preview | Original |
|---|---|---|
| Flower plot | [PNG](flower_plot.png) | [PDF](flower_plot.pdf) |
| Pan/core accumulation | [PNG](pan_core_plot.png) | [PDF](pan_core_plot.pdf) |
| Gene-content MDS | [PNG](gene_content_mds.png) | [PDF](gene_content_mds.pdf) |
| Integrated pan-genome view | [PNG](circos_pangenome.png) | [PDF](circos_pangenome.pdf) |
| BLASTn ANIb | [PNG](ani_blastn_heatmap_strain_names.png) | [PDF](ani_blastn_heatmap_strain_names.pdf) |
| BLASTP POCP | [PNG](pocp_blastp_heatmap_strain_names.png) | [PDF](pocp_blastp_heatmap_strain_names.pdf) |
| Cluster categories per genome | [PNG](category_bar.png) | [PDF](category_bar.pdf) |
| Presence/absence heatmap | [PNG](presence_absence_heatmap.png) | [PDF](presence_absence_heatmap.pdf) |
| Genome Jaccard heatmap | [PNG](genome_jaccard_heatmap.png) | [PDF](genome_jaccard_heatmap.pdf) |
| Gene-cluster intersections | [PNG](euler_upset.png) | [PDF](euler_upset.pdf) |
| Cluster size distribution | [PNG](cluster_size_dist.png) | [PDF](cluster_size_dist.pdf) |
| Sequence identity distribution | [PNG](identity_distribution.png) | [PDF](identity_distribution.pdf) |
| Alignment coverage | [PNG](coverage_scatter.png) | [PDF](coverage_scatter.pdf) |
| Representative protein lengths | [PNG](centroid_length_distribution.png) | [PDF](centroid_length_distribution.pdf) |
| Core-gene conservation | [PNG](core_gene_conservation.png) | [PDF](core_gene_conservation.pdf) |
| Cumulative pan-genome contribution | [PNG](cumulative_pan_contribution.png) | [PDF](cumulative_pan_contribution.pdf) |
| Strain-specific gene contributions | [PNG](singleton_top_n.png) | [PDF](singleton_top_n.pdf) |
| Functional category distribution | [PNG](functional_bar.png) | [PDF](functional_bar.pdf) |
