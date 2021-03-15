library(reticulate)
np <- import("numpy")

corr_mat_file <- [../matrix_files/bw_corr.wg.pearson.npz]

# using reticulate, numpy to import "npz" output from deeptools correlation function
corr_mat_npz <- np$load(corr_mat_file)
corr_mat <- data.frame(corr_mat_npz$f[["matrix"]])

# adding labels and changing order to match the PCC matrix
labels <- gsub("log2_normed_", "", corr_mat_npz$f[["labels"]])
labels <- gsub("_10000.bw", "", labels)
colnames(corr_mat) <- labels
col_order <- c("Ring1b", "H3K27me3", "H2AK119Ub",
              "DNMT3A1_R318W", "DNMT3A1_D333N", "DNMT3A1_K299I", "DNMT3A1_delPWWP-QM", "DNMT3A1_W330R",
              "H3K36me3", "DNMT3A", "H3K36me2")
corr_mat <- corr_mat[, col_order]

##### Plotting

png(filename = "wg_corr_scatterplots.png", height = 2000, width = 2000)
par(mgp = c(3,2,0), cex.axis = 1, cex.lab = 1, col.axis = "white")
pairs(corr_mat, panel=function(x,y){smoothScatter(x,y,add=T)}, lower.panel = NULL, cex.labels = 2, labels=NULL)
dev.off()
