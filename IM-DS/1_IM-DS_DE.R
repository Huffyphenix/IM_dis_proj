library(magrittr)

# path --------------------------------------------------------------------
data_path <- c("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression")
out_path <- c("/project/huff/huff/TKI/result/mRNA_DE/IM_DS")
# load data ---------------------------------------------------------------

all_exp <- readr::read_tsv(file.path(data_path,"AllSamples.GeneExpression.FPKM.xls")) %>%
  dplyr::select(Symbol,G1_1_FPKM,G1_2_FPKM,G1_3_FPKM,G2_1_FPKM,G2_2_FPKM,G2_3_FPKM)
all_exp %>%
  readr::write_rds(file.path(out_path,"all_gene_exp.rds.gz"),compress = "gz")

G2_1 <- readr::read_tsv(file.path(data_path,"G2_1.gene.fpkm.xls")) %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G2_1"="FPKM") %>%
  dplyr::mutate(G2_1=ifelse(G2_1==0,0.01,G2_1))

G2_2 <- readr::read_tsv(file.path(data_path,"G2_2.gene.fpkm.xls")) %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G2_2"="FPKM") %>%
  dplyr::mutate(G2_2=ifelse(G2_2==0,0.01,G2_2))

G2_3 <- readr::read_tsv(file.path(data_path,"G2_3.gene.fpkm.xls")) %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G2_3"="FPKM") %>%
  dplyr::mutate(G2_3=ifelse(G2_3==0,0.01,G2_3))

G1_1 <- readr::read_tsv(file.path(data_path,"G1_1.gene.fpkm.xls")) %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G1_1"="FPKM") %>%
  dplyr::mutate(G1_1=ifelse(G1_1==0,0.01,G1_1))

G1_2 <- readr::read_tsv(file.path(data_path,"G1_2.gene.fpkm.xls")) %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G1_2"="FPKM") %>%
  dplyr::mutate(G1_2=ifelse(G1_2==0,0.01,G1_2))

G1_3 <- readr::read_tsv(file.path(data_path,"G1_3.gene.fpkm.xls")) %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G1_3"="FPKM") %>%
  dplyr::mutate(G1_3=ifelse(G1_3==0,0.01,G1_3))

DE_all <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/DifferentiallyExpressedGene/DEGList/G1-VS-G2.NOIseq_Method.GeneDiffExp.xls")

# 3 sample De overlap -----
fc_threshold <- 0.585

all_exp %>%
  dplyr::mutate(G2_G1_1=ifelse(log2(G2_1_FPKM/G1_1_FPKM)>=fc_threshold,"Up","None")) %>%
  dplyr::mutate(G2_G1_1=ifelse(log2(G2_1_FPKM/G1_1_FPKM)<=(-fc_threshold),"Down",G2_G1_1)) %>%
  dplyr::mutate(G2_G1_2=ifelse(log2(G2_2_FPKM/G1_2_FPKM)>=fc_threshold,"Up","None")) %>%
  dplyr::mutate(G2_G1_2=ifelse(log2(G2_2_FPKM/G1_2_FPKM)<=(-fc_threshold),"Down",G2_G1_2)) %>%
  dplyr::mutate(G2_G1_3=ifelse(log2(G2_3_FPKM/G1_3_FPKM)>=fc_threshold,"Up","None")) %>%
  dplyr::mutate(G2_G1_3=ifelse(log2(G2_3_FPKM/G1_3_FPKM)<=(-fc_threshold),"Down",G2_G1_3)) -> all_samples_DE_info

# get overlap of 3 samples
all_samples_DE_info %>%
  dplyr::filter(G2_G1_1=="Up" & G2_G1_2=="Up" & G2_G1_3=="Up") -> Up_in_all_samples
all_samples_DE_info %>%
  dplyr::filter(G2_G1_1=="Down" & G2_G1_2=="Down" & G2_G1_3=="Down") -> Down_in_all_samples

fn_test <- function(all_samples_DE_info,trend){
  all_samples_DE_info %>%
    dplyr::filter(G2_G1_1==trend & G2_G1_2==trend) -> in_1_2
  all_samples_DE_info %>%
    dplyr::filter(G2_G1_1==trend & G2_G1_3==trend) -> in_1_3
  all_samples_DE_info %>%
    dplyr::filter(G2_G1_2==trend & G2_G1_3==trend) -> in_2_3
  rbind(in_1_2,in_1_3) %>% 
    rbind(in_2_3)  %>%
    unique()
}
all_samples_DE_info %>% fn_test(trend = "Up") -> Up_in_2_samples
all_samples_DE_info %>% fn_test(trend = "Down") -> Down_in_2_samples

# fn_at_least_2<- function(a,b,c){
#   a %>%
#     intersect(b) -> d_1
#   a %>%
#     intersect(c) -> d_2
#   b %>%
#     intersect(c) -> d_3
#   c(d_1,d_2,d_3) %>% unique() ->result
#   return(result)
# }
# fn_at_least_2(DE_1_up$Symbol,DE_2_up$Symbol,DE_3_up$Symbol) -> genes_up_in_at_least_2_pairs
# fn_at_least_2(DE_1_down$Symbol,DE_2_down$Symbol,DE_3_down$Symbol) -> genes_down_in_at_least_2_pairs
# 
# c(genes_up_in_all_pairs,genes_down_in_all_pairs)-> genes_DE_in_at_least_2_pairs

# noiseq result overlap ---------------------------------------------------
# BGI result filter
# filter condition: FC1.5 probability0.8
DE_all %>%
  dplyr::select(Symbol,`G2-Expression`,`G1-Expression`,`log2FoldChange(G2/G1)`,Probability) %>%
  dplyr::filter(Probability>=0.8) %>%
  dplyr::filter(abs(`log2FoldChange(G2/G1)`)>=fc_threshold) %>%
  dplyr::rename("log2FC"=`log2FoldChange(G2/G1)`) %>%
  dplyr::mutate(`DS/IM`=ifelse(log2FC>0,"Up","Down")) -> BGI_DE_all_0.8_1.5

BGI_DE_all_0.8_1.5 %>%
  dplyr::inner_join(all_samples_DE_info,by="Symbol") -> BGI_DE_all_0.8_1.5.info
BGI_DE_all_0.8_1.5.info %>%
  readr::write_tsv(path = file.path(out_path,"DS-IM_BGI_DE_0.8_1.5_mRNA.info"))


BGI_DE_all_0.8_1.5 %>%
  dplyr::select(Symbol) %>%
  dplyr::inner_join(all_exp,by="Symbol") -> BGI_DE_all_0.8_1.5_exp
BGI_DE_all_0.8_1.5_exp %>%
  readr::write_tsv(path = file.path(out_path,"BGI_DE_all_0.8_1.5.exp"))
# overlap with 3 sample
# all samples
BGI_DE_all_0.8_1.5 %>%
  dplyr::inner_join(Down_in_all_samples,by="Symbol") -> Down_in_all_test
BGI_DE_all_0.8_1.5 %>%
  dplyr::inner_join(Up_in_all_samples,by="Symbol") -> Up_in_all_test
rbind(Down_in_all_test,Up_in_all_test) %>%
  readr::write_tsv(path = file.path(out_path,"DE_in_all-AND-in_BGI_0.8_1.5"))
c(Down_in_all_test$Symbol,Up_in_all_test$Symbol) -> DE_in_all_test.list

# at least 2 samples
BGI_DE_all_0.8_1.5 %>%
  dplyr::inner_join(Down_in_2_samples,by="Symbol") -> Down_in_2_and_all_test
BGI_DE_all_0.8_1.5 %>%
  dplyr::inner_join(Up_in_2_samples,by="Symbol") -> Up_in_2_and_all_test
rbind(Down_in_2_and_all_test,Up_in_2_and_all_test) %>%
  readr::write_tsv(path = file.path(out_path,"DE_in_2-AND-in_BGI_0.8_1.5"))
c(Down_in_2_and_all_test$Symbol,Up_in_2_and_all_test$Symbol) -> DE_in_2_and_all_test.list

# statistic --------------------------------------------------------------
BGI_DE_all_0.8_1.5$`DS/IM` %>% table()


# plot --------------------------------------------------------------------

DE_all%>%
  dplyr::rename("log2FC"=`log2FoldChange(G2/G1)`) %>%
  dplyr::mutate(`G2/G1`=ifelse(log2FC>=(0.585) ,"Up","None")) %>%
  dplyr::mutate(`G2/G1`=ifelse(log2FC<=(-0.585),"Down",`G2/G1`)) %>%
  dplyr::mutate(`G2/G1`=ifelse(Probability>=0.8,`G2/G1`,"None")) %>%
  dplyr::mutate(`G2/G1`=ifelse(is.na(Probability),"None",`G2/G1`)) %>%
  dplyr::select(Symbol,log2FC,Probability,`G2/G1`) %>%
  dplyr::mutate(alpha=ifelse(Symbol %in% DE_in_2_and_all_test.list & `G2/G1` !="None",0.5,0.1)) %>%
  dplyr::mutate(alpha=ifelse(Symbol %in% DE_in_all_test.list & `G2/G1` !="None",1,alpha)) %>%
  # dplyr::mutate(Probability=-log10(1-Probability)) %>%
  dplyr::mutate(color=ifelse(`G2/G1`=="Up","red","grey")) %>%
  dplyr::mutate(color=ifelse(`G2/G1`=="Down","blue",color)) -> point_ready

library(ggplot2)

# FC and Probability distribution
point_ready %>%
  ggplot() +
  geom_point(aes(x=log2FC,y=Probability,color=`G2/G1`,alpha=alpha)) + #,colour=point_ready$color
  scale_alpha_continuous(
    name="Significant Group",
    limits=c(0.1,1),
    breaks=c(0.1,0.5,1),
    labels=c("Only Noiseq","At 2 samples & Noiseq","All 3 samples & Noiseq")
  ) +
  xlab("Log2(FC)") +
  theme(
    legend.position = 'bottom',
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 17),
    axis.text = element_text(size = 17),
    legend.text = element_text(size = 17),
    legend.title = element_text(size = 20)
  )-> p;p

ggsave(file.path(out_path,"plot/G1-control.MA.plot.pdf"),p,device = "pdf",width = 3,height = 3)

# Heatmap of DE genes ---------
library(ComplexHeatmap)

# exp data prepare
DE_exp <- as.data.frame(BGI_DE_all_0.8_1.5_exp %>% dplyr::filter(Symbol %in% DE_in_2_and_all_test.list))
rownames(DE_exp) <- DE_exp$Symbol  
DE_exp <- DE_exp[,-1] %>% as.matrix()
colnames(DE_exp) <- sub("_FPKM","",colnames(DE_exp))

# annotation
BGI_DE_all_0.8_1.5 %>%
 dplyr::filter(Symbol %in% DE_in_2_and_all_test.list) -> DE.info
DE.info %>%
  dplyr::select(`DS/IM`) %>%
  as.data.frame() -> DE_anno
rownames(DE_anno) <- DE.info$Symbol  

# row annotation
gene_anno_plot = rowAnnotation(df = DE_anno,
                    col = list(`DS/IM` = c("Up" = "red", "Down" = "green")),
                    width = unit(0.5, "cm")
)

sam_anno <- data.frame(Group = c(rep("IM",3),rep("DS",3)))
rownames(sam_anno) <- DE_exp %>% colnames()
DE_exp[DE_exp==0] <- 0.01
log2(DE_exp)->log2DE_exp
sam_anno_plot = HeatmapAnnotation(df = sam_anno,
                        boxplot = anno_boxplot(log2DE_exp, axis = TRUE),
                        col = list(Group = c("DS" = "pink", "IM" = "purple")))


DE_exp %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() -> DE_exp_rowscale
pdf(file.path(out_path,"plot/DE_mRNA_exp_heatmap.pdf"),width = 5,height = 6)
pdf(file.path(out_path,"plot/DE_all_test_mRNA_exp_heatmap.pdf"),width = 5,height = 6)
pdf(file.path(out_path,"plot/DE_2_sample_and_noiseq_mRNA_exp_heatmap.pdf"),width = 5,height = 6)
he = Heatmap(DE_exp_rowscale, 
             show_row_names = TRUE, 
             cluster_columns = FALSE,
             top_annotation = sam_anno_plot, top_annotation_height = unit(3, "cm"),
             heatmap_legend_param = list(title = c("Experssion")))
he + gene_anno_plot
dev.off()
