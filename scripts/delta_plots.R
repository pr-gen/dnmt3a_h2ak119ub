library(rtracklayer)
library(ggplot2)
library(svglite)

##### Import delta data
h2ak119ub_delta <- import(con = "[../data/delta/log2_normed_H2AK119Ub_10000_delta-subtract_forced.bedgraph]",
                          format = "BED")
DNMT3A1_delPWWP_QM_delta <- import(con = "[../data/delta/log2_normed_DNMT3A1_delPWWP-QM_10000_delta-subtract_forced.bedgraph",
                                   format = "BED")
DNMT3A1_K299I_delta <- import(con = "[../data/delta/log2_normed_DNMT3A1_K299I_10000_delta-subtract_forced.bedgraph",
                              format = "BED")
DNMT3A1_W330R_delta <- import(con = "[../data/delta/log2_normed_DNMT3A1_W330R_10000_delta-subtract_forced.bedgraph",
                              format = "BED")
DNMT3A1_R318W_delta <- import(con = "[../data/delta/log2_normed_DNMT3A1_R318W_10000_delta-subtract_forced.bedgraph",
                              format = "BED")
DNMT3A1_D333N_delta <- import(con = "[../data/delta/log2_normed_DNMT3A1_D333N_10000_delta-subtract_forced.bedgraph",
                              format = "BED")

delta_h2ak119ub <- as.numeric(h2ak119ub_delta$name)  
delta_DNMT3A1_delPWWP_QM <- as.numeric(DNMT3A1_delPWWP_QM_delta$name)
delta_DNMT3A1_K299I <- as.numeric(DNMT3A1_K299I_delta$name)
delta_DNMT3A1_W330R <- as.numeric(DNMT3A1_W330R_delta$name)
delta_DNMT3A1_R318W <- as.numeric(DNMT3A1_R318W_delta$name)
delta_DNMT3A1_D333N <- as.numeric(DNMT3A1_D333N_delta$name)

df <- data.frame(delta_h2ak119ub, delta_DNMT3A1_delPWWP_QM, delta_DNMT3A1_K299I, delta_DNMT3A1_W330R, delta_DNMT3A1_R318W, delta_DNMT3A1_D333N)


smooth_plot <- function(h2ak119ub, DNMT3A1_mut, mut_name) {
  outdir <- "[../plot_files]"
  pcc_h2ak119ub <- cor(x = h2ak119ub, y = DNMT3A1_mut, method = "pearson", use = "complete.obs")
  
  h2ak119ub_lab <- "delta H2AK119Ub\n(sgRing1a/b vs parental)"
  DNMT3A1_mut_lab <- paste("delta DNMT3A1 ", mut_name, "\n(sgRing1a/b vs parental)", sep = "")
  
  svg(filename = paste(mut_name,"delta_h2ak119ub_genome-wide_scatterplot.svg",sep="_")) 
  mar.default <- c(5,4,4,2) + 0.1
  par(mgp=c(3.5,1,0))
  par(mar = mar.default + c(1,2, 0, 0))
  smoothScatter(x = h2ak119ub, y = DNMT3A1_mut, 
                nbin = 128,
                colramp = colorRampPalette(c('white', '#1a5d8f')),
                nrpoints = 100, ret.selection = FALSE,
                pch = ".", cex = 1, col = "black",
                postPlotHook = box,
                xlab = h2ak119ub_lab, ylab = DNMT3A1_mut_lab,
                xaxs = par("xaxs"), yaxs = par("yaxs"),); legend("topright", paste("PCC: ", signif(pcc_h2ak119ub, digits=3), sep=""), bty="n")
  dev.off()
  
}

smooth_plot(df$delta_h2ak119ub, df$delta_h3k27me3, df$delta_DNMT3A1_delPWWP_QM, "delPWWP")
smooth_plot(df$delta_h2ak119ub, df$delta_h3k27me3, df$delta_DNMT3A1_K299I, "K299I")
smooth_plot(df$delta_h2ak119ub, df$delta_h3k27me3, df$delta_DNMT3A1_W330R, "W330R")
smooth_plot(df$delta_h2ak119ub, df$delta_h3k27me3, df$delta_DNMT3A1_R318W, "R318W")
smooth_plot(df$delta_h2ak119ub, df$delta_h3k27me3, df$delta_DNMT3A1_D333N, "D333N")
```
