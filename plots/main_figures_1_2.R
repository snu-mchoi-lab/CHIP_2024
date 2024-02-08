# updated in 2024.02.08. 
# Hyungtai Sim rpendity@outlook.com


# required library
vlib = c("tidyverse", "data.table", "ggsci", "ggpubr", "lemon", "Matrix")
lapply(vlib, require, character.only = TRUE, quietly = TRUE) |> suppressMessages()

## check working directories
## setwd("/data/podo/Projects/project_HS/202111-CHIP/CHIPv4_batch5")

# 1. preparing sample information

## the information contains some private information. to retrieve this, contact us.
cleaned_si_long = read_delim("230906_joint_clinical_info_v5_cleaned_long.txt", sep = "\t")

cleaned_si_long = cleaned_si_long %>% 
  mutate(
  # Create Bins
  Age_Bin = floor(Age/5)*5, 
  category_VAF = ifelse(CHIP_MaxAF >= 0.1, "VAF >= 0.1",
                  ifelse(CHIP_MaxAF > 0.05, "0.1 > VAF >= 0.05", 
                    ifelse(CHIP_MaxAF >= 0.02, "0.05 > VAF >= 0.02", "Negative"))),
  ) %>% 
  mutate(
  # Create Levels
  Pathology = factor(Pathology, levels = c("LUAD", "LUSC", "OTHERS")),
  Response = factor(Response, levels = c("PR", "PD", "SD")),
  time = factor(time, levels = c("before", "after")),
  category_VAF = factor(category_VAF, levels = c("Negative", "0.05 > VAF >= 0.02", "0.1 > VAF >= 0.05", "VAF >= 0.1"))
  )

# 2. import variants
df_variants = fread("230328_annotation_merged_CHIP.txt", sep = "\t") %>%
  left_join(., cleaned_si_long, by = c("ID"="CHIP_ID")) %>%
  filter(is.na(bindID)==FALSE) %>% 
  mutate(tracing_ID = paste(bindID, CHROM, POS, REF, ALT, Gene_Name, sep = "-"))  %>% 
  select(1:5,9:12,Annotation, Gene_Name, HGVS.p, HGVS.c,Existing_variation, IO_ID:tracing_ID) %>%
  group_by(tracing_ID) %>% mutate(var_n = n()) 

df_variants_raw = fread("230324_annotation_merged_raw_output.txt", sep = "\t") %>%
  left_join(., cleaned_si_long, by = c("ID"="CHIP_ID")) %>%
  filter(is.na(bindID)==FALSE) %>% 
  mutate(tracing_ID = paste(bindID, CHROM, POS, REF, ALT, Gene_Name, sep = "-"))  %>% 
  select(1:5,9:12,Annotation, Gene_Name, HGVS.p, HGVS.c,Existing_variation, IO_ID:tracing_ID) %>%
  group_by(tracing_ID) 

df_variants_enlarged = df_variants_raw %>% filter(tracing_ID %in% df_variants$tracing_ID) %>%
  filter(gt_AF >= 0.01) %>%
  mutate(var_n = n())

# 3. F1b

p1b = df_variants %>% 
  mutate(time = factor(time, levels = c("before", "after"))) %>% 
  group_by(Gene_Name) %>% mutate(gene_n = n()) %>% arrange(-gene_n) %>% filter(gene_n >=2) %>%
  mutate(Gene_Name = factor(Gene_Name, levels = unique(.$Gene_Name))) %>%
  ggplot(aes(x = Gene_Name, fill = time)) + geom_bar(position = position_dodge2(preserve = "single")) + 
  scale_y_continuous(breaks= c(2,4,6,8,10))+
  scale_fill_manual(values = c("#3F72AF", "#B83B5E")) +
  theme_pubr() + rotate_x_text(angle = 45)+
  xlab("Gene Name") + ylab("Count")

ggsave("../figure_prep/p1b.pdf", p1b, width = 2200, height = 1200, units = "px", dpi = 300)


# 4. F1c

# stats

df_variants_enlarged %>%
  filter(var_n == 2) %>%
  mutate(time = factor(time, levels = c("before", "after"))) %>%
  select(tracing_ID, gt_AF, time, Pathology, Gene_Name) %>%
  pivot_wider(names_from = time, values_from = gt_AF) %>% lm(before ~ after, data  = .) %>% 
  summary() %>% 
  broom::tidy()

p1c = df_variants_enlarged %>%
  filter(var_n == 2) %>%
  mutate(time = factor(time, levels = c("before", "after"))) %>%
  select(tracing_ID, gt_AF, time, Pathology, Gene_Name) %>%
  pivot_wider(names_from = time, values_from = gt_AF) %>%
  ggplot(aes(x = before, y = after)) + 
    geom_point(aes(color = Pathology), size = 3) +
    geom_smooth(method = "lm", se = F, linetype = "dashed", color = "grey40", alpha = 0.5)+
    theme_pubr() +
    xlab("Allele Frequency, Before Treatment") +
    ylab("Allele Frequency, Before Treatment") +
    scale_color_manual(values = c("#F3AA60","#468B97","grey40")) +
    stat_cor(label.x = 0.12, label.y = 0.2) +
    scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25), limits = c(0, 0.27))+
    scale_y_continuous(breaks=  c(0, 0.05, 0.10, 0.15, 0.20), limits =c(0, 0.23))

ggsave("../figure_prep/p1c.pdf", p1b, width = 1500, height = 1500, units = "px", dpi = 300)

# 5. F1d

p1d1 = df_variants_enlarged %>%
  filter(var_n == 2) %>%
  mutate(time = factor(time, levels = c("before", "after")), name = "All genes") %>%
  ggplot(aes(x = time, y = gt_AF)) + 
  geom_line(aes(group = tracing_ID, color = Pathology), alpha = 0.8, size = 0.90) +
  geom_point(aes(fill = Age), shape = 21, size = 3, stroke= 1.5)+
  theme_pubr(legend = "none")+
  scale_color_manual(values = c("#004C99","#00994C","grey20")) +
  scale_fill_gradient(low = "#008dD6", high= "grey70")+
  scale_y_log10(breaks= c(0, 0.1, 0.5), limits = c(0.01, 0.5))+
  annotation_logticks(sides = "l")+
  facet_wrap(~name)+
  theme(strip.background = element_rect(fill= "#FFFFFF", color = "#FFFFFF"),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 13))+
  xlab("") +
  ylab("VAF of  CH variants")

p1d2 = df_variants_enlarged %>%
  filter(var_n == 2, Gene_Name %in% c("DNMT3A", "PPM1D", "TET2")) %>%
  mutate(time = factor(time, levels = c("before", "after"))) %>%
  ggplot(aes(x = time, y = gt_AF)) + 
  geom_line(aes(group = tracing_ID, color = Pathology), alpha = 0.8, size = 0.90) +
  geom_point(aes(fill = Age), shape = 21, size = 3, stroke= 1.5)+
  theme_pubr(legend = "right")+
  scale_color_manual(values = c("#004C99","#00994C","grey20")) +
  scale_fill_gradient(low = "#008dD6", high= "grey70")+
  scale_y_log10(breaks= c(0, 0.1, 0.5), limits = c(0.01, 0.5))+
  annotation_logticks(sides = "l")+
  xlab("Sampling Time") +
  ylab("")+
  theme(panel.border=element_blank(), 
        axis.line=element_line(linewidth=0.5),
        strip.background = element_rect(fill= "#FFFFFF", color = "#FFFFFF"),
        strip.text = element_text(size = 13, face = "italic"),
        axis.title.x = element_text(hjust = 0.29),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_blank()
        )+
  facet_wrap(~Gene_Name, scale = 'free')

ggsave("../figure_prep/p1d.pdf", 
       p1d1 + p1d2 + plot_layout(widths = c(1, 3.1)), 
       width = 9, height = 6, units = "in", dpi = 300)


# 6. F2a

p2a = cleaned_si_long %>%
  filter(Pathology %in% c("LUAD", "LUSC")) %>%
  group_by(Pathology, category_VAF) %>%
  summarise(n = n()) %>%
  mutate(pos_y = rev(cumsum(rev(n)/sum(n)))) %>%
  mutate(box = "ALL") %>%
  rename("VAF Bins" = category_VAF) %>%
  ggplot(aes(x = Pathology, y = n, fill =  `VAF Bins`)) + 
  geom_bar(stat = "identity", position  = "fill") +
  scale_fill_manual(values = rev(c("#bb2200","#bb6633","#ddaa88","grey70")))+
  geom_text(aes(label = n, y = pos_y), vjust = 1.5) +
  facet_wrap(~box)+
  theme_pubr(legend = "none") +
  xlab("Pathology") +
  ylab("N")

ggsave("../figure_prep/p2a.pdf", p2a, width = 800, height = 1800, units = "px", dpi = 400)

# 7. F2b

p2b = cleaned_si_long %>%
  filter(Pathology %in% c("LUAD", "LUSC")) %>%
  group_by(Pathology,time, category_VAF) %>%
  summarise(n = n()) %>%
  mutate(pos_y = rev(cumsum(rev(n))/sum(n))) %>%
  rename("VAF Bins" = category_VAF) %>%
  ggplot(aes(x = time, y = n, fill =  `VAF Bins`))+ 
  geom_bar(stat = "identity", position  = "fill") +
  scale_fill_manual(values = rev(c("#bb2200","#bb6633","#ddaa88","grey70")))+
  geom_text(aes(label = n, y = pos_y), vjust = 1.5) +
  facet_wrap(~Pathology)+
  theme_pubr(legend = "right") +
  xlab("Time") +
  ylab("N")

ggsave("../figure_prep/p2b.pdf", p2b, width = 2000, height = 1800, units = "px", dpi = 400)

# 8. F2c

p2c = cleaned_si_long %>%
  filter(Pathology %in% c("LUAD", "LUSC"),
         Response %in% c("PR", "PD")) %>%
  group_by(Pathology,Response, category_VAF) %>%
  summarise(n = n()) %>%
  mutate(pos_y = rev(cumsum(rev(n))/sum(n))) %>%
  rename("VAF Bins" = category_VAF) %>%
  ggplot(aes(x = Response, y = n, fill =  `VAF Bins`)) + 
  geom_bar(stat = "identity", position  = "fill") +
  scale_fill_manual(values = rev(c("#bb2200","#bb6633","#ddaa88","grey70")))+
  geom_text(aes(label = n, y = pos_y), vjust = 1.5) +
  facet_wrap(~Pathology)+
  theme_pubr(legend = "right") +
  xlab("Response") +
  ylab("N")

ggsave("../figure_prep/p2c.pdf", p2c, width = 2000, height = 1800, units = "px", dpi = 400)

# 9. F2d

p2d = cleaned_si_long %>%
  filter(Pathology %in% c("LUSC", "LUAD"),
         Response %in% c("PR", "PD")) %>%
  group_by(Pathology, Age_Bin, category_VAF) %>%
  summarise(n = n()) %>%
  mutate(pos_y = rev(cumsum(rev(n))/sum(n))) %>%
  rename("VAF Bins" = category_VAF) %>%
  ggplot(aes(x = Age_Bin, y = n, fill =  `VAF Bins`)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Pathology, ncol = 4) +
  scale_fill_manual(values = rev(c("#bb2200","#bb6633","#ddaa88","grey70")))+
  geom_text(aes(label = n, y = pos_y), vjust = 1.5) +
  theme_pubr(legend = "right") +
  xlab("Age") +
  ylab("N")

ggsave("../figure_prep/p2d.pdf", p2b1, width = 2400, height = 1800, units = "px", dpi = 400)


## WES related figures

chr_anno_path = "/data/podo/Projects/project_HS/202111-CHIP/CHIP_PBMCWES/annotation/"

chr_raw_file = "../230710_annotation_merged_raw_output.txt"
chr_filtered_file = "../230919_annotation_merged_putativeCHIP_v2.txt"
chr_filtered_file = "../230919_annotation_merged_panelmatched.txt"

df_filtered_file = fread(chr_filtered_file)%>% 
  group_by(bind) %>%
  add_count() %>%
  ungroup() %>%
  select(ID, CHROM, POS, gt_AF, gt_DP, 
         Annotation, Gene_Name, HGVS.c, HGVS.p,Existing_variation, n)

gene_set = fread("../../CHIP/SGI_CHIP/Gene_List_SGI_CHIP_V1_Rev2_TE-93684413_hg19.txt", header = F)
gene_set_lCHIP = fread("../../CHIP_TcWES_batch2/230405_gene_list_lchip.txt", header = F) 

gene_set = rbind(gene_set, gene_set_lCHIP) %>% distinct()

df_si = fread("../230712_sample_information.txt")

df_variants = df_filtered_file %>% group_by(ID) %>%
  filter(Gene_Name %in% gene_set$V1, gt_AF < 0.35 | gt_AF > 0.65, gt_DP >= 40) %>%
  summarise(CHIP_N = n(), 
            CHIP_MaxAF = max(gt_AF),
            CHIP_binary = "TRUE",
            CHIP_Severity = ifelse(max(gt_AF) >= 0.1 | n() >= 2, "TRUE", "FALSE"))

df_anno_variants = left_join(df_si, df_variants, by = "ID") %>%
  mutate(CHIP_N = ifelse(is.na(CHIP_N) == FALSE, CHIP_N, 0),
         CHIP_MaxAF = ifelse(is.na(CHIP_MaxAF) == FALSE, CHIP_MaxAF, 0),
         CHIP_binary = ifelse(is.na(CHIP_binary) == FALSE, CHIP_binary, FALSE),
         CHIP_Severity = ifelse(is.na(CHIP_Severity) == FALSE, CHIP_Severity, FALSE)) %>%
  mutate(
    # Create Bins
    Age_Bin = floor(Age/5)*5, 
    category_VAF = ifelse(CHIP_MaxAF >= 0.1, "VAF >= 0.1",
                          ifelse(CHIP_MaxAF > 0.05, "0.1 > VAF >= 0.05", 
                                 ifelse(CHIP_MaxAF >= 0.02, "0.05 > VAF >= 0.02", "Negative"))),
  ) %>% mutate(
    # Create Levels
    Pathology = factor(Histology, levels = c("LUAD", "LUSC", "OTHERS")),
    Response = factor(Response, levels = c("PR", "PD", "SD")),
    category_VAF = factor(category_VAF, levels = c("Negative", "0.05 > VAF >= 0.02", "0.1 > VAF >= 0.05", "VAF >= 0.1"))
  ) %>% select(-Histology)
df_anno_variants

# p2e
p2e = df_anno_variants %>%
  filter(Pathology %in% c("LUAD", "LUSC")) %>%
  group_by(Pathology, category_VAF) %>%
  summarise(n = n()) %>%
  mutate(pos_y = rev(cumsum(rev(n))/sum(n)),
         NAME = "ALL") %>%
  rename("VAF Bins" = category_VAF) %>%
  ggplot(aes(x = Pathology, y = n, fill =  `VAF Bins`)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = rev(c("#bb2200","#bb6633","#ddaa88","grey70")))+
  geom_text(aes(label = n, y = pos_y), vjust = 1.5) +
  facet_wrap(~NAME)+
  theme_pubr(legend = "none") +
  xlab("Pathology") +
  ylab("Proportion")

ggsave("figure_prep/p2e.pdf", p2e, width = 800, height = 1800, units = "px", dpi = 400)

# p2f

p2f = df_anno_variants %>%
  filter(Pathology %in% c("LUAD", "LUSC"),
         Response %in% c("PR", "PD")) %>%
  group_by(Pathology,Response, category_VAF) %>%
  summarise(n = n()) %>%
  mutate(pos_y = rev(cumsum(rev(n))/sum(n))) %>%
  rename("VAF Bins" = category_VAF) %>%
  ggplot(aes(x = Response, y = n, fill =  `VAF Bins`)) + 
  geom_bar(stat = "identity", position  = "fill") +
  scale_fill_manual(values = rev(c("#bb2200","#bb6633","#ddaa88","grey70")))+
  geom_text(aes(label = n, y = pos_y), vjust = 1.5) +
  facet_wrap(~Pathology)+
  theme_pubr(legend = "right") +
  xlab("Response") +
  ylab("N")

ggsave("figure_prep/p2f.pdf", p2c2_f, width = 2000, height = 1800, units = "px", dpi = 400)

# p2g

p2g = df_anno_variants %>%
  mutate(Age_bin = floor(Age/5)*5) %>%
  filter(Pathology %in% c("LUAD", "LUSC")) %>%
  group_by(Age_bin, Pathology, category_VAF) %>%
  summarise(n = n()) %>%
  mutate(pos_y = rev(cumsum(rev(n))/sum(n)),
         NAME = "ALL") %>%
  rename("VAF Bins" = category_VAF) %>%
  ggplot(aes(x = Age_bin, y = n, fill = `VAF Bins`)) + 
  geom_bar(stat = "identity", position  = "fill") +
  scale_fill_manual(values = rev(c("#bb2200","#bb6633","#ddaa88","grey70")))+
  scale_x_continuous(breaks = seq(30, 80, by = 10))+
  geom_text(aes(label = n, y = pos_y), vjust = 1.5) +
  facet_wrap(~Pathology)+
  theme_pubr(legend = "none") +
  xlab("Age") + ylab("")

p2g
ggsave("figure_prep/p2g.pdf", p2g, width = 2300, height = 1800, units = "px", dpi = 400)
