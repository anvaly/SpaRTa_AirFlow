# Author: Anna Lyubetskaya. Date: 21-03-05
# Following: https://irrationone.github.io/cellassign/articles/introduction-to-cellassign.html


install.packages("tensorflow")
library(tensorflow)


## ERROR HERE ----


install_tensorflow(extra_packages = "tensorflow-probability")


## DIDN'T TEST ----


tensorflow::tf_config()

sess = tf$Session()
hello <- tf$constant('Hello, TensorFlow!')
sess$run(hello)

BiocManager::install('cellassign')

devtools::install_github("Irrationone/cellassign")

library(SingleCellExperiment)
library(cellassign)

data(example_sce)
print(example_sce)

print(head(example_sce$Group))

data(example_marker_mat)
print(example_marker_mat)

s <- sizeFactors(example_sce)

fit <- cellassign(exprs_obj = example_sce[rownames(example_marker_mat),], 
                  marker_gene_info = example_marker_mat, 
                  s = s, 
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)

print(fit)

print(head(celltypes(fit)))

print(str(mleparams(fit)))

pheatmap::pheatmap(cellprobs(fit))

print(table(example_sce$Group, celltypes(fit)))
