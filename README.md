# Trying to Replicate GSEA: Subramanian et al. (2005)

This project is trying to replicate the Gene Set Enrichment Analysis (GSEA) methodology as described in the landmark article:

- Subramanian, A., Tamayo, P., Mootha, V.K., et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43), 15545-15550. [https://www.pnas.org/doi/10.1073/pnas.0506580102](https://www.pnas.org/doi/10.1073/pnas.0506580102)

## Data Sources

All data used in this project are publicly available from the Broad Institute's GSEA/MSigDB resource:
- GSEA Example Datasets: [https://www.gsea-msigdb.org/gsea/datasets.jsp](https://www.gsea-msigdb.org/gsea/datasets.jsp)

Specifically, for the leukemia analysis:
- **Expression data:** [Leukemia_hgu95av2.gct](https://data.broadinstitute.org/gsea-msigdb/gsea/dataset_files/Leukemia_hgu95av2.gct)
- **Gene sets:** [c2.symbols.gmt](https://data.broadinstitute.org/gsea-msigdb/gsea/dataset_files/c2.symbols.gmt)

## Workflow

1. **Compute a ranked gene list** using the signal-to-noise ratio between ALL and AML samples from the leukemia dataset.
2. **Run a from-scratch GSEA pre-ranked analysis** using the ranked list and the C2 (curated) gene sets.
3. **Compare results** to those reported in the original GSEA article.

## Current Issue

**All tested gene sets return ES=0 and p=1.0.** This likely indicates a mismatch between gene symbols in the ranked list and those in the gene sets, or another data integration issue. Resolving this is a current priority for the project.

## Files
- `compute_snr_ranked_list.py`: Computes the signal-to-noise ratio ranked gene list from the .gct and .cls files.
- `gsea_preranked_with_fdr_and_plot.py`: Runs the GSEA pre-ranked analysis and outputs results and enrichment plots.
- `leukemia_ranked_gene_list.txt`: Ranked gene list for GSEA.
- `gsea_preranked_results_with_fdr.csv`: GSEA results (ES, p-value, FDR/q-value).
- `gsea_plots/`: Enrichment plots for selected gene sets.

## References
- [Original GSEA article (PNAS, 2005)](https://www.pnas.org/doi/10.1073/pnas.0506580102)
- [GSEA Example Datasets](https://www.gsea-msigdb.org/gsea/datasets.jsp)
- [Leukemia_hgu95av2.gct](https://data.broadinstitute.org/gsea-msigdb/gsea/dataset_files/Leukemia_hgu95av2.gct)
- [c2.symbols.gmt](https://data.broadinstitute.org/gsea-msigdb/gsea/dataset_files/c2.symbols.gmt)
