library(annotatr)
library(rtracklayer)
library(dplyr)
library(ggplot2)

##### Import log2 normalized ChIP-seq data for H2AK119Ub
# This script describes H2AK119Ub (Ext. Data Fig. 5e). The corresponding figure for  H3K27me3 (Ext. Data Fig. 5f)
# was produced identically by swapping files for normalized H3K27me3 ChIP-seq bedgraphs and changing Rx values.
# Ext. Data Fig. 2d was produced using this script as well (without Rx normalization).

h2ak119ub_par <- read_regions(con = "[../data/log2_normed_parental/log2_normed_H2AK119Ub_10000_forced.bedgraph",
                              genome = 'mm10', rename_name = "score", format = "BED")
h2ak119ub_sgEzh2 <- read_regions(con = "[../data/log2_normed_sgEzh2/log2_normed_H2AK119Ub_10000_forced.bedgraph]",
                                 genome = 'mm10', rename_name = "score", format = "BED")
h2ak119ub_sgRing1 <- read_regions(con = "[../data/log2_normed_sgRing1/log2_normed_H2AK119Ub_10000_forced.bedgraph]",
                                  genome = 'mm10', rename_name = "score", format = "BED")

##### Select annotations for intersection with regions

read_annotations(con = "[../data/H3K36me2_peaks.bed]",
                 genome = "mm10", name = "h3k36me2_peaks", format = "BED")
read_annotations(con = "[../data/H3K27me3_peaks.bed]",
                 genome = "mm10", name = "h3k27me3_peaks", format = "BED")
read_annotations(con = "[../data/consensus_H2AK119Ub_peaks.sorted.bed]",
                 genome = "mm10", name = "h2ak119ub_sicer_con_peaks", format = "BED")

# Inclusion/exclusion regions

read_annotations(con = "[../data/H2AUb-peaks_excluding_H3K27me3-peaks.bed]",
                 genome = "mm10", name = "h2ak119ub_no_h3k27me3", format = "BED")

read_annotations(con = "[../data/H2AUb-peaks_overlapping_H3K27me3-peaks.bed]",
                 genome = "mm10", name = "h2ak119ub_yes_h3k27me3", format = "BED")

read_annotations(con = "[../data/H3K27me3-peaks_excluding_H2AK119Ub-peaks.bed]",
                 genome = "mm10", name = "h3k27me3_no_h2ak119ub", format = "BED")

# RefSeq genes
read_annotations(con = "[../data/mm10_refseq_genes_merged.bed]",
                 genome = "mm10", name = "genebodies", format = "BED")

annots = c(
  'mm10_cpgs', 'mm10_basicgenes', ### built-in
  'mm10_custom_h3k36me2_peaks', 'mm10_custom_h3k27me3_peaks', 
  'mm10_custom_h2ak119ub_sicer_peaks', 'mm10_custom_h2ak119ub_sicer_con_peaks', 'mm10_custom_h2ak119ub_no_h3k27me3', 'mm10_custom_h2ak119ub_yes_h3k27me3', 'mm10_custom_h3k27me3_no_h2ak119ub',
  'mm10_custom_genebodies') ### custom (see annotatr docs for naming convention)

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'mm10', annotations = annots)

# Intersect the regions we read in with the annotations

h2ak119ub_par_annotated = annotate_regions(
  regions = h2ak119ub_par,
  annotations = annotations, ignore.strand = TRUE, quiet = FALSE)
h2ak119ub_sgEzh2_annotated = annotate_regions(
  regions = h2ak119ub_sgEzh2,
  annotations = annotations, ignore.strand = TRUE, quiet = FALSE)
h2ak119ub_sgRing1_annotated = annotate_regions(
  regions = h2ak119ub_sgRing1,
  annotations = annotations, ignore.strand = TRUE, quiet = FALSE)
# Coerce GRanges objects into dfs
df_h2ak119ub_par_annotated = data.frame(h2ak119ub_par_annotated)
df_h2ak119ub_par_annotated$name = factor(c("Parental"))
df_h2ak119ub_par_annotated$score <- as.numeric(2^as.numeric(df_h2ak119ub_par_annotated$score)) # converting log2 ratio to ratio
df_h2ak119ub_sgEzh2_annotated = data.frame(h2ak119ub_sgEzh2_annotated)
df_h2ak119ub_sgEzh2_annotated$name = factor(c("sgEzh2"))
df_h2ak119ub_sgEzh2_annotated$score <- as.numeric(2^as.numeric(df_h2ak119ub_sgEzh2_annotated$score)) # converting log2 ratio to ratio
df_h2ak119ub_sgRing1_annotated = data.frame(h2ak119ub_sgRing1_annotated)
df_h2ak119ub_sgRing1_annotated$name = factor(c("sgRing1a/b"))
df_h2ak119ub_sgRing1_annotated$score <- as.numeric(2^as.numeric(df_h2ak119ub_sgRing1_annotated$score)) # converting log2 ratio to ratio

### Rx factor normalization relative to parental ###
rx_sgEzh2_par <- 1.119/1.243 # sgEzh2 vs par Rx ratio
df_h2ak119ub_sgEzh2_annotated$score <- as.numeric(df_h2ak119ub_sgEzh2_annotated$score) * rx_sgEzh2_par

rx_sgRing1_par <- 0.1961/1.243 # sgRing1 vs par Rx ratio
df_h2ak119ub_sgRing1_annotated$score <- as.numeric(df_h2ak119ub_sgRing1_annotated$score) * rx_sgRing1_par


df <- rbind(df_h2ak119ub_par_annotated, df_h2ak119ub_sgRing1_annotated, df_h2ak119ub_sgEzh2_annotated)

##### Summarizing the data
slimmed_df <- data.frame(df$score, df$annot.type, df$name)
slimmed_df$df.score <- as.numeric(as.character(slimmed_df$df.score))
slimmed_df <- slimmed_df[complete.cases(slimmed_df), ]
sampled_slimmed_df <- slimmed_df[sample(nrow(slimmed_df), 1000), ]
summary_stats <- slimmed_df %>%
  group_by(df.annot.type, df.name) %>% 
  dplyr::summarize(avg = mean(df.score))

print(tbl_df(summary_stats), n=100)

##### Extracting annotations of interest and specifying order
annots_order <- c('mm10_custom_genebodies', 'mm10_genes_promoters', 'mm10_cpg_islands',
                  'mm10_custom_h2ak119ub_sicer_con_peaks', 
                  'mm10_custom_h3k36me2_peaks', 'mm10_custom_h3k27me3_peaks',
                  'mm10_custom_h2ak119ub_no_h3k27me3', 'mm10_custom_h2ak119ub_yes_h3k27me3', 'mm10_custom_h3k27me3_no_h2ak119ub')
annots_labels <- c('Gene bodies', 'Promoters', 'CGIs', 'H2AK119Ub\npeaks', 'H3K36me2\npeaks', 'H3K27me3\npeaks', 'H2AK119Ub,\nK27me3-' ,'H2AK119Ub,\nK27me3+','H3K27me3,\nH2AK119Ub-')

annots_stats <- summary_stats %>%
  filter(df.annot.type %in% annots_order)

# Setting order of annotations
annots_stats$df.annot.type <- factor(annots_stats$df.annot.type, levels = rev(annots_order))

##### Plotting
p <- ggplot(annots_stats, aes(x = df.annot.type, y = avg, fill = df.name)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.75) +
  coord_flip() +
  xlab(element_blank()) +
  scale_x_discrete(breaks = annots_order, labels = annots_labels) +
  ylab("Ratio(observed/expected)") +
  ggtitle("H2AK119Ub")

p <- p +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 1, color = "black") +
  theme(panel.grid.major.x = element_line(colour = "grey")) + 
  theme(panel.grid.major.y = element_blank()) +
  scale_fill_manual(values=c("#999999", "#cb0000", "#800080")) +
  theme(legend.position="bottom")
print(p)