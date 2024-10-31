vlib = c("tidyverse", "data.table", "Seurat", "tidyseurat", "SeuratDisk",
         "ggsci", "ggpubr", "WGCNA", "hdWGCNA", "patchwork",
         "GeneOverlap", "enrichR", "scales", "fgsea", "msigdbr")

# need seurat 4 for executing following codes.
# this object contains data from wgcna analysis.

seurat_obj = readRDS("230808_hdWGCNA.RDS")

update_md = fread("230809_all_merged.txt.gz")
colnames(update_md) = c("cell", "rn_anno_l1", "rn_anno_l2")

temp_md = seurat_obj@meta.data %>% rownames_to_column("cell") 

new_md = left_join(temp_md, update_md) %>%
  mutate(anno_l1 = ifelse(is.na(rn_anno_l1) ==FALSE, rn_anno_l1, anno_l1)) %>%
  mutate(anno_l2 = ifelse(is.na(rn_anno_l2) ==FALSE, rn_anno_l2, anno_l2)) %>%
  select(-rn_anno_l1, -rn_anno_l2) %>%
  mutate(VAF_bin = ifelse(CHIP_MaxAF >= 0.1, "VAF>=0.1",
                          ifelse(CHIP_MaxAF >= 0.02, "0.1>VAF>=0.02", "Negative"))) %>%
  mutate(VAF_bin= factor(VAF_bin, levels = c("VAF>=0.1", "0.1>VAF>=0.02", "Negative"))) %>%
  mutate(time = factor(time, levels = c("base", "1st"))) %>%
  mutate(group = paste(sep = "-", VAF_bin, anno_l2, time, Pathology)) %>%
  column_to_rownames("cell")

seurat_obj %>% saveRDS("230808_hdWGCNA_updatemd.RDS")

## DE analysis

anno_l2_all = seurat_obj@meta.data %>% arrange(anno_l2) %>% pull(anno_l2) %>% unique()

list_DE_res_base = vector(mode = "list", length = length(anno_l2_all))
list_DE_res_1st = vector(mode = "list", length = length(anno_l2_all))
names(list_DE_res_base) = anno_l2_all
names(list_DE_res_1st) = anno_l2_all

DE_by_cluster_base = function(seurat_obj, chr_group_by, chr_assay = "RNA", cluster_name){
  seurat_obj_subset = seurat_obj %>% filter(anno_l2 == cluster_name, time == "base")
  DE_results = presto::wilcoxauc(seurat_obj_subset, 
                                 group_by = chr_group_by,
                                 seurat_assay = chr_assay)
  return(DE_results)
}

DE_by_cluster_1st = function(seurat_obj, chr_group_by, chr_assay = "RNA", cluster_name){
  seurat_obj_subset = seurat_obj %>% filter(anno_l2 == cluster_name, time == "1st")
  DE_results = presto::wilcoxauc(seurat_obj_subset, 
                                 group_by = chr_group_by,
                                 seurat_assay = chr_assay)
  return(DE_results)
}

for (each_idx in 1:length(anno_l2_all)){
  each_cluster = anno_l2_all[each_idx]
  each_DE_results_base = DE_by_cluster_base(seurat_obj, chr_group_by = "group", chr_assay = "SCT", cluster_name = each_cluster)
  each_DE_results_1st = DE_by_cluster_1st(seurat_obj, chr_group_by = "group", chr_assay = "SCT", cluster_name = each_cluster)
  list_DE_res_base[[each_idx]] = each_DE_results_base
  list_DE_res_1st[[each_idx]] = each_DE_results_1st
}

list_DE_res_base %>% saveRDS("CHIP_scRNA/230803_DE_list_base_RNA.RDS")
list_DE_res_1st %>% saveRDS("CHIP_scRNA/230803_DE_list_1st_RNA.RDS")

for (each_idx in 1:length(anno_l2_all)){
  list_DE_res[[each_idx]] = list_DE_res[[each_idx]] %>% filter(padj < 1e-5 & (auc < 0.45 | auc > 0.55))
}

## GSEA analysis

m_df = msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets = m_df %>% split(x = .$gene_symbol, f = .$gs_name)


chr_names = names(list_DE_res_base)[1:25]
list_res_gsea_base = vector(mode = "list", length = 5L)
names(list_res_gsea_base) = c("VAF>=0.1-LUSC", "0.1>VAF>=0.02-LUSC", "Negative-LUSC", "0.1>VAF>=0.02-LUAD", "Negative-LUAD")
list_res_gsea_base$`VAF>=0.1-LUSC` = vector(mode = "list", length = 25L)
list_res_gsea_base$`0.1>VAF>=0.02-LUSC` = vector(mode = "list", length = 25L)
list_res_gsea_base$`Negative-LUSC` = vector(mode = "list", length = 25L)
list_res_gsea_base$`0.1>VAF>=0.02-LUAD` = vector(mode = "list", length = 25L)
list_res_gsea_base$`Negative-LUAD` = vector(mode = "list", length = 25L)

names(list_res_gsea_base$`VAF>=0.1-LUSC`) = chr_names 
names(list_res_gsea_base$`0.1>VAF>=0.02-LUSC`) = chr_names 
names(list_res_gsea_base$`Negative-LUSC`) = chr_names 
names(list_res_gsea_base$`0.1>VAF>=0.02-LUAD`) = chr_names 
names(list_res_gsea_base$`Negative-LUAD`) = chr_names 

list_rank_base = vector(mode = "list", length = 5L)
names(list_rank_base) = c("VAF>=0.1-LUSC", "0.1>VAF>=0.02-LUSC", "Negative-LUSC", "0.1>VAF>=0.02-LUAD", "Negative-LUAD")
list_rank_base$`VAF>=0.1-LUSC` = vector(mode = "list", length = 25L)
list_rank_base$`0.1>VAF>=0.02-LUSC` = vector(mode = "list", length = 25L)
list_rank_base$`Negative-LUSC` = vector(mode = "list", length = 25L)
list_rank_base$`0.1>VAF>=0.02-LUAD` = vector(mode = "list", length = 25L)
list_rank_base$`Negative-LUAD` = vector(mode = "list", length = 25L)
names(list_rank_base$`VAF>=0.1-LUSC`)= chr_names 
names(list_rank_base$`0.1>VAF>=0.02-LUSC`)=chr_names 
names(list_rank_base$`Negative-LUSC`)= chr_names 
names(list_rank_base$`0.1>VAF>=0.02-LUAD`)= chr_names 
names(list_rank_base$`Negative-LUAD`)= chr_names 

for (each_name in chr_names) {
  print(paste0(each_name, " cluster is processing"))
  list_gsea = list_DE_res_base[[each_name]] %>% 
    separate(group, sep = "-", into = c("VAF_bin", "cluster", "time", "pathology")) %>%
    mutate(group = paste0(VAF_bin,"-", pathology)) %>%
    filter(pathology %in% c("LUAD", "LUSC"), avgExpr > 0.3) %>%
    mutate(group = factor(group, levels = c("VAF>=0.1-LUSC", "0.1>VAF>=0.02-LUSC", "Negative-LUSC", "0.1>VAF>=0.02-LUAD", "Negative-LUAD"))) %>%
    group_by(group) %>% 
    arrange(-auc, -abs(logFC)) %>%
    select(group, feature, auc) 
  
  c_names_list_gsea = list_gsea %>% group_keys() %>% pull(group)
  list_gsea = list_gsea %>% group_split(.keep = F)
  names(list_gsea) = c_names_list_gsea
  
  
  df_gsea_stat = list_DE_res_base[[each_name]] %>% 
    separate(group, sep = "-", into = c("VAF_bin", "cluster", "time", "pathology")) %>%
    mutate(group = paste0(VAF_bin,"-", pathology)) %>%
    filter(pathology %in% c("LUAD", "LUSC")) %>%
    group_by(group) %>% 
    arrange(-auc, -abs(logFC)) %>%
    select(group, feature, logFC, auc, padj)
  
  p1 = df_gsea_stat %>% 
    mutate(group = factor(group, levels = c("VAF>=0.1-LUSC", "0.1>VAF>=0.02-LUSC", "Negative-LUSC", "0.1>VAF>=0.02-LUAD", "Negative-LUAD"))) %>%
    ggplot(aes(x = auc, fill = group)) + 
    geom_histogram(binwidth = 0.005) + 
    scale_fill_igv() +
    scale_y_log10()+
    theme_pubr(legend = "none") +
    xlim(0.3,0.8)
  
  p2 = df_gsea_stat %>% 
    mutate(group = factor(group, levels = c("VAF>=0.1-LUSC", "0.1>VAF>=0.02-LUSC", "Negative-LUSC", "0.1>VAF>=0.02-LUAD", "Negative-LUAD"))) %>%
    filter(padj < 1e-3) %>%
    ggplot(aes(x = logFC, y = -log10(padj))) +
    scale_color_igv() +
    geom_point(aes(color = group)) +
    theme_pubr(legend = "bottom")
  
  
  for (idx in 1:5){
    #rank =  list_gsea[[idx]] %>% filter(!auc == 0.5) %>% mutate(auc = auc - 0.5) %>% deframe()
    rank =  list_gsea[[idx]] %>% filter(auc > 0.6 | auc < 0.4) %>% mutate(auc = auc - 0.5) %>% deframe()
    fgsea_res = fgsea(fgsea_sets,rank, eps = 0)
    list_rank_base[[idx]][[each_name]] = rank
    list_res_gsea_base[[idx]][[each_name]] = fgsea_res
  }
}