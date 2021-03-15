library(ggplot2)
library(gridExtra)
library(rtracklayer)
library(reshape2)
library(ggpubr)
library(plyr)

##### Import data #####
# Each data point for wt_in_cgi and mut_in_cgi represents the ChIP-seq signal of a 10-kb bin overlapping CpG islands  
# Each data point for wt_in_random and mut_in_random is the the mean ChIP-seq signal across a random 10-kb bin within the genome
# "overlapping" files generated using UCSC CpG islands, and bedtools overlap as described in Methods
wt_in_cgi <- import(con = "[../data/log2_normed_parental/log2_normed_DNMT3A1_10000_forced_overlapping_CGIs_bl_removed.bed]",
                    format = "BED")
wt_in_random <- import(con = "[../data/log2_normed_parental/log2_normed_DNMT3A1_10000_forced_mean_per_CGIs_bl_removed_shuffled.bed]",
                       format = "BED")
mut_in_cgi <- import(con = "[../data/log2_normed_parental/log2_normed_DNMT3A1_delPWWP_10000_forced_overlapping_CGIs_bl_removed.bed]",
                    format = "BED")
mut_in_random <- import(con = "[../data/log2_normed_parental/log2_normed_DNMT3A1_delPWWP_10000_forced_mean_per_CGIs_bl_removed_shuffled.bed]",
                       format = "BED")

# Isolate scores from  bedgraphs, and form into a dataframe
wt_in_cgi <- as.numeric(wt_in_cgi$name)
wt_in_random <- as.numeric(wt_in_random$name)
mut_in_cgi <- as.numeric(mut_in_cgi$name)
mut_in_random <- as.numeric(mut_in_random$name)
df <- data.frame(wt_in_cgi, wt_in_random, mut_in_cgi, mut_in_random)

##### Format data
meth_df <- data.frame(melt(na_dropped_df))
meth_df$value <- meth_df$value 
meth_df$region <- NA
meth_df$mutation <- NA
meth_df[meth_df$variable == 'wt_in_cgi' | meth_df$variable == 'mut_in_cgi',]$region <- c("CGI")
meth_df[meth_df$variable == 'wt_in_random' | meth_df$variable == 'mut_in_random',]$region <- c("Random")
meth_df[meth_df$variable == 'wt_in_cgi' | meth_df$variable == 'wt_in_random',]$mutation <- c("WT")
meth_df[meth_df$variable == 'mut_in_cgi' | meth_df$variable == 'mut_in_random',]$mutation <- c("delPWWP")
meth_df$mutation <- factor(meth_df$mutation, levels = c("WT", "delPWWP"))
meth_df$region <- factor(meth_df$region, levels = c("Random",  "CGI"))
# -----------

##### Density plots 
# These are designed based on the ones from Weinberg et al. 2019
# Building a single density plot
p1 <- ggplot(meth_df, aes(x=value, fill=region)) + 
  geom_density(color="transparent", alpha=0.5)
p1 <- p1 + scale_fill_manual(values=c("#999999", "#E69F00"))
p1 <- p1 + theme_classic()

# Faceting density plots by mutation type
p <- p1 + facet_grid(. ~ mutation) +
  theme(strip.background = element_rect(
    color="transparent", fill="transparent", size=1.5, linetype="solid"
  ))
p <- p +
  xlab("Signal/input ratio (log2)") +
  ylab("Density") +
  labs(fill = "Region")


mu <- ddply(meth_df, c("mutation", "region"), summarise, grp.mean=mean(value))
head(mu)
p <- p + 
  geom_vline(data=mu, aes(xintercept=grp.mean, color=region),
                linetype="dashed", show.legend = FALSE) +
  scale_color_manual(values=c("#999999", "#E69F00")) 

print(p)
