vlib = c("tidyverse", "mashr", "ashr", "ggpubr", "data.table", "patchwork", "cowplot",
         "future.apply", "arrow", 'pheatmap', "Seurat", "hdWGCNA", "magrittr", "optparse")
lapply(vlib, require, character.only = TRUE, quietly = TRUE) |> suppressMessages()

base.dir = "/path/to/base/"

setwd(paste0(base.dir, "202310-sceqtl_v6/"))

saveRDS.gz <- function(object,file,threads=4) {
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}
readRDS.gz <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

set.seed(42)

# run with 
# Rscript path/to/thisfile -c Mono

option_list = list(
  make_option(c("-c", "--cluster"), action="store", default="NK", type='character',
              help="select a cluster for analysis.")
  
)

opt = parse_args(OptionParser(option_list=option_list))
c_name = opt$cluster

# optionally enable multithreading
enableWGCNAThreads(nThreads = 4)
seurat_obj = readRDS.gz("scRNAseq_downstream_scenic/seurat_obj_preped_for_wgcna.RDS")


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c_name, # the name of the group of interest in the group.by column
  group.by='anno_l1', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)


wrap_plots(plot_list, ncol=2) %>% ggsave(paste0("./scRNAseq_downstream_scenic/hdWGCNA_soft_power_", c_name ,".png"), .,
                                         height = 3000, width = 3000, units = "px", dpi = 300)

seurat_obj %>% saveRDS.gz(file = paste0("./scRNAseq_downstream_scenic/seurat_obj_wgcna_", c_name, ".RDS"))


# define soft-power of each clusters

soft = 10 # check the soft power using soft_power png.


seurat_obj = readRDS.gz(file = paste0( "./scRNAseq_downstream_scenic/seurat_obj_wgcna_",c_name,".RDS"))


seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = opt$soft,
  detectCutHeight = 0.980, # 995
  mergeCutHeight = 0.05, # 05
  tom_name = c_name,
  overwrite_tom =TRUE  # name of the topoligical overlap matrix written to disk
)

png(paste0("scRNAseq_downstream_scenic/Dendrogram_", c_name, ".png"), width = 4000, height = 2400, res = 600)
seurat_obj %>% PlotDendrogram()
dev.off()

seurat_obj = seurat_obj %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="sample_time"
)

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'anno_l1', group_name = c_name,
)

seurat_obj %>% saveRDS.gz(file = paste0( "./scRNAseq_downstream_scenic/seurat_obj_wgcna_",c_name,".RDS"))

