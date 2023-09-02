library(tidyverse)

# read in result files
eigenValues <- read_delim("pca3.eigenval", delim = " ", col_names = F)
eigenVectors <- read_delim("pca3.eigenvec", delim = "\t", col_names = T)

## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

# PCA plot
 
ggplot(data = eigenVectors) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = IID, shape = IID), size = 3, show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of selected populations",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"),
       colour = 'IID', shape = 'IID') +
  theme_minimal()


https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/pca/
