
library (smartsnp, lib.loc = "~/_lib")

My <- c(rep("A",2), rep("B",2))
sm.pca <- smart_pca(snp_data = "Z_filtered.traw",
                    missing_value = NA,
                        sample_group = My)
pdf("my_plot.pdf")
plot(sm.pca$pca.sample_coordinates[, c(3,4)])
dev.off()
 
