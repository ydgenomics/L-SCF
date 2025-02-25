# /software/conda/Anaconda/bin/R 250223
# https://github.com/cellgeni/sceasy
library(optparse)
library(sceasy)
library(reticulate)
library(Seurat)
packageVersion("Seurat")
option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = NULL,
    help = "Path to input file for convrting"
  ),
  make_option(c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "Path to output file for convrting"
  ),
  make_option(c("-a", "--assay"),
    type = "character", default = "RNA",
    help = "Assay name for the output file"
  ),
  make_option(c("-m", "--main_layer"),
    type = "character", default = "RNA",
    help = "Main layer name for the output file"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

use_python("/software/conda/Anaconda/bin/python")

loompy <- reticulate::import('loompy')

temp0 <- readRDS(opt$input_file)

print(temp0)
DefaultAssay(temp0) <- 'RNA'
print(temp0)

# Check '_index' in meta.data and meta.feature, because RLIGER way will product '_index' cause error of python dataframe
if ("_index" %in% colnames(temp0@meta.data)) {
  print("Change _index into new_index for meta.data")
  colnames(temp0@meta.data)[colnames(temp0@meta.data) == "_index"] <- "new_index"
}
if ("meta.features" %in% slotNames(temp0[["RNA"]])) {
  features_meta <- temp0[["RNA"]]@meta.features
  if ("_index" %in% colnames(features_meta)) {
    colnames(features_meta)[colnames(features_meta) == "_index"] <- "new_index"
    print("Change _index into new_index for meta.features")
  }
  temp0[["RNA"]]@meta.features <- features_meta
} else {
  message("meta.features does not exist in the RNA Assay of temp0.")
}
#
colnames(temp0@meta.data)
#
temp0[["RNA"]] <- as(temp0[["RNA"]], "Assay")
#
sceasy::convertFormat(temp0, from="seurat", to="anndata", assay = opt$assay, main_layer=opt$main_layer, outFile = opt$output_file)
