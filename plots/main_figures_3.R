vlib = c("tidyverse", "data.table", "Seurat", "tidyseurat", "SeuratDisk",
         "ggsci", "ggpubr", "WGCNA", "hdWGCNA", "patchwork",
         "GeneOverlap", "enrichR", "scales", "fgsea", "msigdbr")
lapply(vlib, require, character.only = TRUE, quietly = TRUE) |> suppressMessages()
setwd("/data/podo/Projects/project_HS/202111-CHIP/")

c_colors = c(
    "#9cc964","#1077f3", "#5ba8f7", "#e83326", "#f98517", 
    "#fbaa5d", "#c85b00", "#ff7d67", "#d3ba59", "#eee8b6",
    "#997600", "#c6a000", "#175182","#48846c", "#ac0000", 
    "#ee74ee", "#970098", "#cc34cd", "#7018d3", "#9b54f3", 
    "#cab1e5", "#903900", "#0050ae", "#008c5c", "grey50",
    "#83ca9c"
    )

seurat_obj = readRDS("230808_hdWGCNA_updatemd.RDS")


# F3a
p3a = DimPlot(seurat_obj, 
             group.by = "anno_l2", reduction = "umap", cols = alpha(c_colors, 0.66),
             raster.dpi = c(4096,4096), pt.size = 8, label = T) + 
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        plot.title = element_blank(), 
        legend.position =  "none")

# F3b legend

p3b1L = seurat_obj@meta.data %>% filter(Pathology %in% c("LUAD", "LUSC")) %>%
  mutate(time = ifelse(time == "base", "before", "after" )) %>%
  mutate(time = factor(time, levels = c("before", "after"))) %>%
  ggplot(aes(x = "1", fill = anno_l2)) + geom_bar(position = "fill", color = "black", alpha = 0.66) + scale_fill_manual(values = c_colors) + 
  xlab(' ')+ylab('proportion')+
  theme_pubr(legend = "right") 

p3b1L =get_legend(p3b1L, position = NULL) %>% as_ggplot() 

ggsave("figure_prep/p3b1L.pdf", p3b1L,  height = 4, width = 3, units = "in", dpi = 300)

# F3b full

p3b1 = seurat_obj@meta.data %>% filter(Pathology %in% c("LUAD", "LUSC")) %>%
  mutate(time = ifelse(time == "base", "before", "after" )) %>%
  mutate(time = factor(time, levels = c("before", "after"))) %>%
  ggplot(aes(x = time, fill = anno_l2)) + geom_bar(position = "fill", color = "black", alpha = 0.66) + scale_fill_manual(values = c_colors) + 
  xlab(' ')+ylab('proportion')+
  theme_pubr(legend = "none") + facet_wrap(~Pathology)

p3b2 = seurat_obj@meta.data %>% filter(Pathology %in% c("LUAD", "LUSC")) %>%
  mutate(time = ifelse(time == "base", "before", "after" )) %>%
  mutate(time = factor(time, levels = c("before", "after"))) %>%
  mutate(VAF_bin = factor(VAF_bin, levels = c("Negative","0.1>VAF>=0.02",  "VAF>=0.1"))) %>%
  ggplot(aes(x = time, fill = anno_l2)) + geom_bar(position = "fill", color = "black", alpha = 0.66) + scale_fill_manual(values =  c_colors) + theme_pubr(legend = "none") + facet_wrap(~VAF_bin) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())

gc()
ggsave("figure_prep/p3b1.pdf", p3b1 + p3b2 + plot_layout(widths = c(1,1.5)),  height = 4, width = 5, units = "in", dpi = 300)


list_res_gsea_base = readRDS("./CHIP_scRNA/GSEA/230811_gsea_before_treatment.RDS")
list_res_gsea_1st = readRDS("./CHIP_scRNA/GSEA/230811_gsea_after_treatment.RDS")
list_rank_base = readRDS("./CHIP_scRNA/GSEA/230811_DE_rank_before_treatment.RDS")
list_rank_1st = readRDS("./CHIP_scRNA/GSEA/230811_DE_rank_after_treatment.RDS")

c_pathway = c("TNFA SIGNALING VIA NFKB","INTERFERON ALPHA RESPONSE",
            "OXIDATIVE PHOSPHORYLATION", "HYPOXIA", "UV RESPONSE UP",
             "EPITHELIAL MESENCHYMAL TRANSITION", "P53 PATHWAY")

c_Myeloid = c("CD14_Mono", "CD16_Mono", "cDC")

# F3C was collected from p3c1-p3c4

list_res_gsea_base$`VAF>=0.1-LUSC` %>% bind_rows() %>% arrange(padj)

p3c1 = list_res_gsea_base$`VAF>=0.1-LUSC` %>% bind_rows(.id = 'cluster') %>% 
  arrange(cluster) %>%
  mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
  arrange(pathway) %>%
  mutate(pathway = str_replace(pathway, "^HALLMARK_", ""),
         pathway = str_replace_all(pathway, "_", " ")) %>%
  filter(pathway %in% c_pathway, cluster %in% c_Myeloid) %>%
  mutate(pathway = factor(pathway, levels = rev(c_pathway)),
         cluster = factor(cluster, levels = c_Myeloid)) %>%
  ggplot(aes(x = cluster, y = pathway)) + 
  geom_point(aes(fill = NES, size = -log10(padj)), shape = 21) +
  scale_fill_gradient2(low = 'grey40', high = '#0459b0', oob=squish, limits = c(-2,2))+
  scale_size_continuous(limits = function(x){c(1, max(x))}, range = c(0,8))+
  theme_pubr(legend = "none", x.text.angle = 45)+
  theme(axis.text.x = element_text(vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave("./figure_prep/p3c1.pdf", p3c1,  width = 1600, height = 1000, units = "px", dpi = 300)

list_res_gsea_base$`Negative-LUSC` %>% bind_rows() %>% arrange(padj)
p3c2 = list_res_gsea_base$`Negative-LUSC` %>% bind_rows() %>% 
  arrange(cluster) %>%
  mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
  arrange(pathway) %>%
  mutate(pathway = str_replace(pathway, "^HALLMARK_", ""),
         pathway = str_replace_all(pathway, "_", " ")) %>%
  filter(pathway %in% c_pathway, cluster %in% c_Myeloid) %>%
  mutate(pathway = factor(pathway, levels = rev(c_pathway)),
         cluster = factor(cluster, levels = c_Myeloid)) %>%
  ggplot(aes(x = cluster, y = pathway)) + 
  geom_point(aes(fill = NES, size = -log10(padj)), shape = 21) +
  scale_fill_gradient2(low = 'grey40', high = '#0459b0', oob=squish, limits = c(-2,2))+
  scale_size_continuous(limits = function(x){c(1, 15)}, range = c(0,8))+
  theme_pubr(legend = "none", x.text.angle = 45)+
  theme(axis.text.x = element_text(vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave("./figure_prep/p3c2.pdf", p3c2,  width = 1600, height = 1000, units = "px", dpi = 300)


p3c3 = list_res_gsea_1st$`VAF>=0.1-LUSC` %>% bind_rows() %>% 
  arrange(cluster) %>%
  mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
  arrange(pathway) %>%
  mutate(pathway = str_replace(pathway, "^HALLMARK_", ""),
         pathway = str_replace_all(pathway, "_", " ")) %>%
  filter(pathway %in% c_pathway, cluster %in% c_Myeloid) %>%
  mutate(pathway = factor(pathway, levels = rev(c_pathway)),
         cluster = factor(cluster, levels = c_Myeloid)) %>%
  ggplot(aes(x = cluster, y = pathway)) + 
  geom_point(aes(fill = NES, size = -log10(padj)), shape = 21) +
  scale_fill_gradient2(low = 'grey40', high = '#0459b0', oob=squish, limits = c(-2,2))+
  scale_size_continuous(limits = function(x){c(1, 15)}, range = c(0,8))+
  theme_pubr(legend = "none", x.text.angle = 45)+
  theme(axis.text.x = element_text(vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave("./figure_prep/p3c3.pdf", p3c3, width = 1600, height = 1000, units = "px", dpi = 300)

p3c4 = list_res_gsea_1st$`Negative-LUSC` %>% bind_rows() %>% 
  arrange(cluster) %>%
  mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
  arrange(pathway) %>%
  mutate(pathway = str_replace(pathway, "^HALLMARK_", ""),
         pathway = str_replace_all(pathway, "_", " ")) %>%
  filter(pathway %in% c_pathway, cluster %in% c_Myeloid) %>%
  mutate(pathway = factor(pathway, levels = rev(c_pathway)),
         cluster = factor(cluster, levels = c_Myeloid)) %>%
  ggplot(aes(x = cluster, y = pathway)) + 
  geom_point(aes(fill = NES, size = -log10(padj)), shape = 21) +
  scale_fill_gradient2(low = 'grey40', high = '#0459b0', oob=squish, limits = c(-2,2))+
  scale_size_continuous(limits = function(x){c(1, 15)}, range = c(0,8))+
  theme_pubr(legend = "none", x.text.angle = 45)+
  theme(axis.text.x = element_text(vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
p3c4
ggsave("./figure_prep/p3c4.pdf", p3c4 , width = 1600, height = 1000, units = "px", dpi = 300)


# f3d

p3d = seurat_obj %>%
  filter(anno_l2 %in% c("CD14_Mono", "CD16_Mono", "cDC"),
        VAF_bin %in% c("VAF>=0.1", "Negative")) %>%
  mutate(group = paste(time, VAF_bin, sep = "-")) %>%
  mutate(group = factor(group, levels=c("base-Negative", "base-VAF>=0.1", "1st-Negative", "1st-VAF>=0.1",))) %>%
  VlnPlot(features = c("IL1B", "FOSB", "NFKBIA", "JUNB", "SOCS3", "NAMPT", "LAPTM5", "ZFP36"), 
          split.by = "group", assay = "RNA", pt.size = 0, cols = c( "grey60", "#EE5533", "grey35", "#BB3322"), ncol = 2) &
  theme(axis.text.x = element_text(angle = 0, hjust= +0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

ggsave("./figure_prep/p3d.pdf", p3d, width = 3000, height = 2600, units = "px", dpi = 300)''


# f3e

## p3e1
chr_pathway = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
chr_cluster = "CD14_Mono"
chr_condition = "Negative-LUSC"

each_res = list_res_gsea_base[[chr_condition]][[chr_cluster]] %>%
  filter(pathway ==  chr_pathway)
each_res$padj

p3raw =plotEnrichment(fgsea_sets[[chr_pathway]], list_rank_base[[chr_condition]][[chr_cluster]])
df_extract = p3raw$data
df_p3data = data.frame(
  Ranked_list_metric =  list_rank_base[[chr_condition]][[chr_cluster]] 
  ) %>%
  rownames_to_column("gene_name") %>%
  mutate(rank = c(1:nrow(.)))

p3e_1 = df_extract %>%
  ggplot(aes(x = x, y = y)) + 
    geom_line(color ="#3377ff") + 
    geom_hline(yintercept = max(df_p1$y), linetype = "longdash", color = "#dd2200")+
    geom_text(label = paste0("p: ",format(each_res$padj,digits = 3)), 
              x = 0.9*max(df_p1$x), y = 0.9*max(df_p1$y), hjust = 1, vjust = 1) +
    geom_text(label = paste0("NES: ",format(each_res$NES,digits = 3)), 
              x = 0.9*max(df_p1$x), y = 0.9*max(df_p1$y), hjust = 1, vjust = 3) +
    ylab("NES")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  ggtitle(chr_pathway)
p3e_2 = df_extract %>% mutate(sum = cumsum(x)/sum(x)) %>%
  ggplot(aes(x = x, y = y))+
  geom_hline(yintercept = 0)+
  geom_segment(aes(x = x, xend = x, y = 1, yend = -0.8))+
  scale_y_continuous()+ theme_pubr()+ ylab('')+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.5,0,0.5), "cm"))
p3e_3 = df_p3data %>% ggplot(aes(x = rank, y = Ranked_list_metric + 0.5)) + 
  geom_bar(stat = "identity", alpha = 0.8)+ theme_pubr() +
  geom_hline(yintercept = 0.5,  linetype = "longdash")+
  ylab("AUC") + 
  theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"))
p3e1 = cowplot::plot_grid(p3e_1, p3e_2, p3e_3, ncol = 1, align = "v", rel_heights = c(5,1,3), greedy = T)
ggsave("figure_prep/p3e1.pdf", p3e1, width = 4, heigh t = 5, dpi = 300)


## p3e2
chr_pathway = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
chr_cluster = "CD14_Mono"
chr_condition = "VAF>=0.1"

each_res = list_res_gsea_base[[chr_condition]][[chr_cluster]] %>%
  filter(pathway ==  chr_pathway)
each_res$padj

p3raw =plotEnrichment(fgsea_sets[[chr_pathway]], list_rank_base[[chr_condition]][[chr_cluster]])
df_extract = p3raw$data
df_p3data = data.frame(
  Ranked_list_metric =  list_rank_base[[chr_condition]][[chr_cluster]] 
  ) %>%
  rownames_to_column("gene_name") %>%
  mutate(rank = c(1:nrow(.)))

p3e_1 = df_extract %>%
  ggplot(aes(x = x, y = y)) + 
    geom_line(color ="#3377ff") + 
    geom_hline(yintercept = max(df_p1$y), linetype = "longdash", color = "#dd2200")+
    geom_text(label = paste0("p: ",format(each_res$padj,digits = 3)), 
              x = 0.9*max(df_p1$x), y = 0.9*max(df_p1$y), hjust = 1, vjust = 1) +
    geom_text(label = paste0("NES: ",format(each_res$NES,digits = 3)), 
              x = 0.9*max(df_p1$x), y = 0.9*max(df_p1$y), hjust = 1, vjust = 3) +
    ylab("NES")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  ggtitle(chr_pathway)
p3e_2 = df_extract %>% mutate(sum = cumsum(x)/sum(x)) %>%
  ggplot(aes(x = x, y = y))+
  geom_hline(yintercept = 0)+
  geom_segment(aes(x = x, xend = x, y = 1, yend = -0.8))+
  scale_y_continuous()+ theme_pubr()+ ylab('')+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.5,0,0.5), "cm"))
p3e_3 = df_p3data %>% ggplot(aes(x = rank, y = Ranked_list_metric + 0.5)) + 
  geom_bar(stat = "identity", alpha = 0.8)+ theme_pubr() +
  geom_hline(yintercept = 0.5,  linetype = "longdash")+
  ylab("AUC") + 
  theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"))
p3e2 = cowplot::plot_grid(p3e_1, p3e_2, p3e_3, ncol = 1, align = "v", rel_heights = c(5,1,3), greedy = T)
ggsave("figure_prep/p3e2.pdf", p3e2, width = 4, height = 5, dpi = 300)

## p3e3
chr_pathway = "HALLMARK_INTERFERON_GAMMA_RESPONSE"
chr_cluster = "cDC"
chr_condition = "Negative-LUSC"

each_res = list_res_gsea_1st[[chr_condition]][[chr_cluster]] %>%
  filter(pathway ==  chr_pathway)
each_res$padj

p3raw =plotEnrichment(fgsea_sets[[chr_pathway]], list_rank_1st[[chr_condition]][[chr_cluster]])
df_extract = p3raw$data
df_p3data = data.frame(
  Ranked_list_metric =  list_rank_1st[[chr_condition]][[chr_cluster]] 
  ) %>%
  rownames_to_column("gene_name") %>%
  mutate(rank = c(1:nrow(.)))

p3e_1 = df_extract %>%
  ggplot(aes(x = x, y = y)) + 
    geom_line(color ="#3377ff") + 
    geom_hline(yintercept = max(df_p1$y), linetype = "longdash", color = "#dd2200")+
    geom_text(label = paste0("p: ",format(each_res$padj,digits = 3)), 
              x = 0.9*max(df_p1$x), y = 0.9*max(df_p1$y), hjust = 1, vjust = 1) +
    geom_text(label = paste0("NES: ",format(each_res$NES,digits = 3)), 
              x = 0.9*max(df_p1$x), y = 0.9*max(df_p1$y), hjust = 1, vjust = 3) +
    ylab("NES")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  ggtitle(chr_pathway)
p3e_2 = df_extract %>% mutate(sum = cumsum(x)/sum(x)) %>%
  ggplot(aes(x = x, y = y))+
  geom_hline(yintercept = 0)+
  geom_segment(aes(x = x, xend = x, y = 1, yend = -0.8))+
  scale_y_continuous()+ theme_pubr()+ ylab('')+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.5,0,0.5), "cm"))
p3e_3 = df_p3data %>% ggplot(aes(x = rank, y = Ranked_list_metric + 0.5)) + 
  geom_bar(stat = "identity", alpha = 0.8)+ theme_pubr() +
  geom_hline(yintercept = 0.5,  linetype = "longdash")+
  ylab("AUC") + 
  theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"))
p3e3 = cowplot::plot_grid(p3e_1, p3e_2, p3e_3, ncol = 1, align = "v", rel_heights = c(5,1,3), greedy = T)
ggsave("figure_prep/p3e3.pdf", p3e3, width = 4, height = 5, dpi = 300)

## p3e4
chr_pathway = "HALLMARK_INTERFERON_GAMMA_RESPONSE"
chr_cluster = "cDC"
chr_condition = "VAF>=0.1"

each_res = list_res_gsea_1st[[chr_condition]][[chr_cluster]] %>%
  filter(pathway ==  chr_pathway)
each_res$padj

p3raw =plotEnrichment(fgsea_sets[[chr_pathway]], list_rank_1st[[chr_condition]][[chr_cluster]])
df_extract = p3raw$data
df_p3data = data.frame(
  Ranked_list_metric =  list_rank_1st[[chr_condition]][[chr_cluster]] 
  ) %>%
  rownames_to_column("gene_name") %>%
  mutate(rank = c(1:nrow(.)))

p3e_1 = df_extract %>%
  ggplot(aes(x = x, y = y)) + 
    geom_line(color ="#3377ff") + 
    geom_hline(yintercept = max(df_p1$y), linetype = "longdash", color = "#dd2200")+
    geom_text(label = paste0("p: ",format(each_res$padj,digits = 3)), 
              x = 0.9*max(df_p1$x), y = 0.9*max(df_p1$y), hjust = 1, vjust = 1) +
    geom_text(label = paste0("NES: ",format(each_res$NES,digits = 3)), 
              x = 0.9*max(df_p1$x), y = 0.9*max(df_p1$y), hjust = 1, vjust = 3) +
    ylab("NES")+
    theme_pubr()+
    theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  ggtitle(chr_pathway)
p3e_2 = df_extract %>% mutate(sum = cumsum(x)/sum(x)) %>%
  ggplot(aes(x = x, y = y))+
  geom_hline(yintercept = 0)+
  geom_segment(aes(x = x, xend = x, y = 1, yend = -0.8))+
  scale_y_continuous()+ theme_pubr()+ ylab('')+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.5,0,0.5), "cm"))
p3e_3 = df_p3data %>% ggplot(aes(x = rank, y = Ranked_list_metric + 0.5)) + 
  geom_bar(stat = "identity", alpha = 0.8)+ theme_pubr() +
  geom_hline(yintercept = 0.5,  linetype = "longdash")+
  ylab("AUC") + 
  theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"))
p3e4 = cowplot::plot_grid(p3e_1, p3e_2, p3e_3, ncol = 1, align = "v", rel_heights = c(5,1,3), greedy = T)
ggsave("figure_prep/p3e4.pdf", p3e4, width = 4, height = 5, dpi = 300)