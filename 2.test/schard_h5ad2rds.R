# convert_format.R 250224
# Â© EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>
# https://github.com/cellgeni/schard
library(optparse)
library(Seurat)
option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = NULL,
    help = "Path to input file for convrting"
  ),
  make_option(c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "Output file after conversion"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if(FALSE){
opt <- list()
opt$input_file <- "/data/input/Files/ResultData/Workflow/W202409020001454/Cer_leaf__soupx_dataget/Cer_leaf__soupx.h5ad"
opt$output_file <- "/data/work/Cer_leaf__soupx.RDS"
opt$stype <- "anndata_to_seurat"
}
library(schard)
input_file <- opt$input_file
output_file <- opt$output_file
stype="anndata_to_seurat"

if(stype == 'anndata_to_seurat'){
    message(paste0("from anndata to seurat, input: ", input_file))
    dt=schard::h5ad2seurat(input_file)
    saveRDS(dt,file=opt$output_file)
} else if (stype == 'seurat_to_anndata'){
    warning("using SeuratDisk and  another script")
}
