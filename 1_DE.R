library(magrittr)

# path --------------------------------------------------------------------

out_path <- c("/project/huff/huff/TKI/result/mRNA_DE/IM_dis")
# load data ---------------------------------------------------------------

all_exp <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/AllSamples.GeneExpression.FPKM.xls") %>%
  dplyr::select(Symbol,Con_1_FPKM,Con_2_FPKM,Con_3_FPKM,G1_1_FPKM,G1_2_FPKM,G1_3_FPKM)
all_exp %>%
  readr::write_rds(file.path(out_path,"all_gene_exp.rds.gz"),compress = "gz")

Con_1 <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/Con_1.gene.fpkm.xls") %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("Con_1"="FPKM")

Con_2 <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/Con_2.gene.fpkm.xls")%>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("Con_2"="FPKM")

Con_3 <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/Con_3.gene.fpkm.xls") %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("Con_3"="FPKM")

G1_1 <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/G1_1.gene.fpkm.xls") %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G1_1"="FPKM")

G1_2 <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/G1_2.gene.fpkm.xls") %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G1_2"="FPKM")

G1_3 <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/G1_3.gene.fpkm.xls") %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G1_3"="FPKM")

DE_all <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/DifferentiallyExpressedGene/DEGList/Control-VS-G1.NOIseq_Method.GeneDiffExp.xls")

# data merge and anlysis ------------------------------------------------------------------
# fileter condition: FC>1.5
# 3 group, overlap
fc_threshold <- 0.585

Con_1 %>%
  dplyr::inner_join(G1_1,by="Symbol") %>%
  dplyr::mutate(FC = log2(G1_1/Con_1)) %>%
  dplyr::filter(FC>=fc_threshold) -> DE_1_up

Con_1 %>%
  dplyr::inner_join(G1_1,by="Symbol") %>%
  dplyr::mutate(FC = log2(G1_1/Con_1)) %>%
  dplyr::filter(FC <= (-fc_threshold)) -> DE_1_down

Con_2 %>%
  dplyr::inner_join(G1_2,by="Symbol") %>%
  dplyr::mutate(FC = log2(G1_2/Con_2)) %>%
  dplyr::filter(FC>=fc_threshold) -> DE_2_up

Con_2 %>%
  dplyr::inner_join(G1_2,by="Symbol") %>%
  dplyr::mutate(FC = log2(G1_2/Con_2)) %>%
  dplyr::filter(FC <= (-fc_threshold)) -> DE_2_down

Con_3 %>%
  dplyr::inner_join(G1_3,by="Symbol") %>%
  dplyr::mutate(FC = log2(G1_3/Con_3)) %>%
  dplyr::filter(FC>=fc_threshold) -> DE_3_up

Con_3 %>%
  dplyr::inner_join(G1_3,by="Symbol") %>%
  dplyr::mutate(FC = log2(G1_3/Con_3)) %>%
  dplyr::filter(FC <= (-fc_threshold)) -> DE_3_down

DE_1_up %>%
  dplyr::select(Symbol) %>%
  dplyr::semi_join(DE_2_up,by="Symbol") %>%
  dplyr::semi_join(DE_3_up,by="Symbol") %>%
  tidyr::drop_na() %>%
  t() %>%
  as.character()-> genes_up_in_all_pairs

DE_1_down%>%
  dplyr::select(Symbol) %>%
  dplyr::semi_join(DE_2_down,by="Symbol") %>%
  dplyr::semi_join(DE_3_down,by="Symbol") %>%
  tidyr::drop_na() %>%
  t() %>%
  as.character()-> genes_down_in_all_pairs

c(genes_up_in_all_pairs,genes_down_in_all_pairs)-> genes_DE_in_all_pairs

# noiseq result overlap ---------------------------------------------------
# BGI result filter
# filter condition: FC1.5 probability0.8
DE_all %>%
  dplyr::select(Symbol,Probability,`log2FoldChange(G1/Control)`) %>%
  dplyr::filter(Probability>=0.8) %>%
  dplyr::filter(abs(`log2FoldChange(G1/Control)`)>=fc_threshold) %>%
  dplyr::rename("log2FC"=`log2FoldChange(G1/Control)`) %>%
  dplyr::mutate(`G1/Control`=ifelse(log2FC>0,"Up","Down")) -> DE_all_0.8_1.5

DE_all_0.8_1.5 %>%
  readr::write_tsv(path = file.path(out_path,"IM_dis_DE_BGI_0.8_1.5_mRNA"))

DE_all_0.8_1.5_overlap_allPair %>%
  readr::write_tsv(path = file.path(out_path,"IM_dis_DE_mRNA"))

# overlap -----------------------------------------------------------------

DE_all_0.8_1.5 %>%
  dplyr::filter(Symbol %in% genes_DE_in_all_pairs) -> DE_all_0.8_1.5_overlap_allPair

DE_all_0.8_1.5_overlap_allPair %>%
  readr::write_tsv(path = file.path(out_path,"IM_dis_DE_mRNA"))


# statistic ---------------------------------------------------------------

DE_all_0.8_1.5_overlap_allPair %>%
  dplyr::select(`G1/Control`)  %>%
  table()

DE_all_0.8_1.5_overlap_allPair %>%
  dplyr::select(Probability,log2FC) %>%
  dplyr::mutate(Probability=log10(1-Probability)) %>%
  

# plot --------------------------------------------------------------------

DE_all%>%
  dplyr::rename("log2FC"=`log2FoldChange(G1/Control)`) %>%
  dplyr::mutate(`G1/Control`=ifelse(log2FC>=(0.585) ,"Up","None")) %>%
  dplyr::mutate(`G1/Control`=ifelse(log2FC<=(-0.585),"Down",`G1/Control`)) %>%
  dplyr::mutate(`G1/Control`=ifelse(Probability>=0.8,`G1/Control`,"None")) %>%
  dplyr::mutate(`G1/Control`=ifelse(is.na(Probability),"None",`G1/Control`)) %>%
  dplyr::select(Symbol,log2FC,Probability,`G1/Control`) %>%
  dplyr::mutate(alpha=ifelse(Symbol %in% genes_DE_in_all_pairs & `G1/Control` !="None",1,0.5)) %>%
  # dplyr::mutate(Probability=-log10(1-Probability)) %>%
  dplyr::mutate(color=ifelse(`G1/Control`=="Up","red","grey")) %>%
  dplyr::mutate(color=ifelse(`G1/Control`=="Down","blue",color)) -> point_ready

# FC and Probability distribution
library(ggplot2)
point_ready %>%
  ggplot() +
  geom_point(aes(x=log2FC,y=Probability,color=`G1/Control`,alpha=alpha)) + #,colour=point_ready$color
  scale_alpha_continuous(
    name="Concensus in 3 sample",
    limits=c(0.5,1),
    breaks=c(0.5,1),
    labels=c("No","Yes")
  ) +
  xlab("Log2(FC)") +
  theme(
    legend.position = 'bottom'
  )-> p;p

ggsave(file.path(out_path,"plot/G1-control.MA.plot.pdf"),p,device = "pdf",width = 3,height = 3)

# Heatmap of DE genes
library(ComplexHeatmap)

# data prepare
DE_all_0.8_1.5_overlap_allPair %>%
  dplyr::select(Symbol) %>%
  dplyr::inner_join(all_exp,by="Symbol") %>%
  as.data.frame()-> DE_exp
rownames(DE_exp) <- DE_exp$Symbol  
DE_exp <- DE_exp[,-1] %>% as.matrix()
colnames(DE_exp) <- sub("_FPKM","",colnames(DE_exp))

# annotation
DE_all_0.8_1.5_overlap_allPair %>%
  dplyr::select(`G1/Control`) %>%
  as.data.frame() -> DE_anno
rownames(DE_anno) <- DE_all_0.8_1.5_overlap_allPair$Symbol  

# row annotation
ha1 = rowAnnotation(df = DE_anno,
                        col = list(`G1/Control` = c("Up" = "red", "Down" = "blue")),
                        width = unit(0.5, "cm")
)


# ha3 = rowAnnotation(link = row_anno_link(at = subset, labels = labels),
#                     col = list(`G1/Control` = c("Up" = "red", "Down" = "blue")),
#                     width = unit(1, "cm")
# )
# column annotation
sam_anno <- data.frame(Group = c(rep("Con",3),rep("IM",3)))
rownames(sam_anno) <- DE_exp %>% colnames()
log2(DE_exp) -> log2DE_exp
ha2 = HeatmapAnnotation(df = sam_anno,
                        boxplot = anno_boxplot(DE_exp, axis = TRUE),
                        # violin = anno_density(DE_exp, type = "violin", 
                        #                       gp = gpar(fill = c("Con" = "pink", "IM" = "purple"))),
                        col = list(Group = c("Con" = "pink", "IM" = "purple")))

# ha3 = HeatmapAnnotation(boxplot = anno_boxplot(log2(DE_exp), axis = TRUE, axis_direction = "reverse"), 
                        # width = unit(2, "cm"))

DE_exp %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() -> DE_exp_rowscale
pdf(file.path(out_path,"plot/DE_mRNA_exp_heatmap.pdf"),width = 5,height = 6)
he = Heatmap(DE_exp_rowscale, 
             show_row_names = FALSE, 
             top_annotation = ha2, top_annotation_height = unit(6, "cm"),
             heatmap_legend_param = list(title = c("Experssion")))
he +ha1
dev.off()
class(exp_heatmap)

# specific gene expression
# ATPase related

DE_exp[c("PAM16","LTF"),] %>%
  Heatmap()
