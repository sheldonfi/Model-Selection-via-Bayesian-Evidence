library(glmnet)
#install.packages('mogavs')
library(mogavs)
data(crimeData)

y <- crimeData$y
X <- crimeData[,-123]

n <- ncol(X)

################################################################################################
### PCA 

## we scale the matrix X so columns have mean=0 and variance=1

scale(X,center = TRUE, scale = TRUE)

## we perform a singular value decomposition of matrix X

svd_X <- svd(X)
V <- svd_X$v #eigenvectors of sample covariance X^T*X/n
svd_X$d

Z <- svd_X$u %*% diag(svd_X$d) #calculate and store the scores

## perform PCA analysis on X

pca_X <- prcomp(X, center=FALSE, scale. =TRUE)
summary(pca_X)

standard_deviations <- pca_X$sdev
prop_var <- (standard_deviations^2) / sum(standard_deviations^2)
prop_var

#############################################################
## scree plot showing the percentage of variance captured by PCs
library(factoextra)

scree_plot <- fviz_eig(pca_X,ncp = 15, addlabels = FALSE, barfill = '#007FFF', barcolor = '#007FFF')
scree_plot +
  ggtitle("Scree Plot") +                # Add a title
  theme_bw() +                           # Use a black-and-white theme
  theme(plot.title = element_text(hjust = 0.5),  # Center the plot title
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        axis.title = element_text(size = 12),   # Increase font size for axis titles
        axis.text = element_text(size = 10),    # Increase font size for axis labels
        panel.grid.major = element_line(color = "lightgray"),  # Add dashed gridlines
        panel.grid.minor = element_blank(),     # Remove minor gridlines
        panel.background = element_rect(fill = "white"))  # Change panel background color

#############################################################
## similarities and dissimilarities plot

fviz_pca_var(pca_X, col.var="contrib", gradient.cols=c("WHITE", "blue", "red"), geom = 'arrow', label = 'none')
