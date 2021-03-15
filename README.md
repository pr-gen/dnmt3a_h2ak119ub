# Delta score calculation and correlation
1. **bigwig_to_bedgraph.sh**: Map scores to 10kb bins.
2. **bedgraph_deltas.sh**: Calculate delta scores between normalized ChIP-seq bedgraphs.
3. **delta_plot.R**: Plot delta values in a scatterplot.

# Genomewide bigwig correlation scatter plots
1. **multi_bigwig_compare.sh**: Using deepTools for genome-wide correlation of normalized ChIP-seq.
2. **bigwig_correlation_scatterplots.R**: Plotting pairwise scatterplots between all samples.

# Enrichment by genomic region
1. **plot_chip_ratio.R**: Plot ratio of real:expected ChIP-seq reads at annotated genomic regions of interest.

# DNMT3A1 enrichment tests
1. **bigwigToBedgraph.sh**:	Map scores to 10kb bins.
3. **shuffle_regions.sh**: Permute CGIs genome-wide to get a set of comparison regions.
4. **permutation_analysis.R**: Plot DNMT3A scores across CGIs and permuted regions.

# H3K36me2-ranked intergenic region plots
1. **bigwig_to_bedgraph.sh**: Map scores to 10kb bins.
2. **score_per_igr.sh**: Compute average score per IGR.
3. **dnmt3a_across_ranked_IGRs.R**: Plot DNMT3A scores across IGRs, ranked by H3K36me2.
