library(ggplot2)
library(gridExtra)
library(rtracklayer)
library(reshape2)
library(ggpubr)

##### Import data #####
# Each data point is the mean log2 ratio (signal/input) across an entire IGR 
wt_in_igr <- import(con = "[../data/log2_normed_parental/bedgraph_files/mean_per_region/log2_normed_DNMT3A_10000_mean_per_igr_bl-removed.bedgraph]",
                    format = "BED")
mut_in_igr <- import(con = "[../data/log2_normed_parental/bedgraph_files/mean_per_region/log2_normed_DNMT3A1_delPWWP_10000_mean_per_igr_bl-removed.bedgraph",
                     format = "BED")
h3k36me2_in_igr <- import(con = "[../data/log2_normed_parental/bedgraph_files/mean_per_region/log2_normed_H3K36me2_10000_mean_per_igr_bl-removed.bedgraph",
                          format = "BED")

# Isolate scores from  bedgraphs, and form into a dataframe
# All data is in the parental background
wt_in_igr <- as.numeric(wt_in_igr$name)
mut_in_igr <- as.numeric(mut_in_igr$name)
h3k36me2_in_igr <- as.numeric(h3k36me2_in_igr$name)
df <- data.frame(h3k36me2_in_igr, wt_in_igr, mut_in_igr)

##### Keep only complete cases (ChIP-seq normalized log2 score is not N/A for each region across all 4 regions)
na_dropped_df <- df[complete.cases(df), ]
melted_df <- data.frame(melt(na_dropped_df))
melted_df$value <- melted_df$value 
melted_df$region <- NA
melted_df$mark <- NA
melted_df[melted_df$variable == 'wt_in_igr' | melted_df$variable == 'mut_in_igr' | melted_df$variable == 'h3k36me2_in_igr',]$region <- c("IGR")
melted_df[melted_df$variable == 'wt_in_igr',]$mark <- c("WT")
melted_df[melted_df$variable == 'mut_in_igr',]$mark <- c("delPWWP")
melted_df[melted_df$variable == 'h3k36me2_in_igr',]$mark <- c("H3K36me2")
melted_df$mark <- factor(melted_df$mark, levels = c("WT", "delPWWP", "H3K36me2"))
# -----------

# sort by H2AK119Ub
df <- df[order(h3k36me2_in_igr),]

# retain only entries with h3k36me2 > 0 (signal greater than input... i.e. somewhat enriched for k36me2)
df <- df[df$h3k36me2_in_igr > 0,]

### wild type
p <- ggplot(df, aes(x=seq_along(h3k36me2_in_igr)), cex.axis = 2) + 
  geom_point(aes(y=wt_in_igr), colour = "#999999", alpha=0.2, size=0.2) +
  geom_smooth(aes(y=wt_in_igr), method = "loess", color="#999999") +
  labs(x = "IGRs ranked by H3K36me2", y = "DNMT3A1") +
  ylim(-1,2) +
  scale_colour_manual(values = c("#999999", "#E69F00"))

p <- p + scale_y_continuous(sec.axis = sec_axis(~(( . + 0.5)*1.5), name = "H3K36me2")) +
  geom_line(aes(y=h3k36me2_in_igr/1.5-0.5), color="black", size=1, linetype = "solid" ) 

p <- p +
  theme_classic() +
  theme(axis.text = element_text(size = 12))
print(p)

### mutant
q <- ggplot(df, aes(x=seq_along(h3k36me2_in_igr)), cex.axis = 2) + 
  geom_point(aes(y=mut_in_igr), colour = "#E69F00", alpha=0.2, size=0.2) +
  geom_smooth(aes(y=mut_in_igr), method = "loess", color="#E69F00") +
  labs(x = "IGRs ranked by H3K36me2", y = "DNMT3A1") +
  ylim(-1,2) +
  scale_colour_manual(values = c("#999999", "#E69F00"))

q <- q + scale_y_continuous(sec.axis = sec_axis(~(( . + 0.5)*1), name = "H3K36me2")) +
  geom_line(aes(y=h3k36me2_in_igr/1-0.5), color="black", size=1, linetype = "solid" ) 

q <- q +
  theme_classic() +
  theme(axis.text = element_text(size = 12))
print(q)

