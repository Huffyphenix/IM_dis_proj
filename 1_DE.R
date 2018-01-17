library(magrittr)

# path --------------------------------------------------------------------

out_path <- c("/project/huff/huff/TKI/result/mRNA_DE")
# load data ---------------------------------------------------------------

all_exp <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/AllSamples.GeneExpression.FPKM.xls") %>%
  dplyr::select(Symbol,Con_1_FPKM,Con_2_FPKM,Con_3_FPKM,G1_1_FPKM,G1_2_FPKM,G1_3_FPKM)
  
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

G1_2 <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/G1_1.gene.fpkm.xls") %>%
  dplyr::select(Symbol,FPKM) %>%
  dplyr::rename("G1_2"="FPKM")

G1_3 <- readr::read_tsv("/project/huff/huff/TKI/data/RNA/F17FTSCCWLJ2064_HUMwcyR/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression/G1_1.gene.fpkm.xls") %>%
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
  t() %>%
  as.character()-> genes_up_in_all_pairs

DE_1_down%>%
  dplyr::select(Symbol) %>%
  dplyr::semi_join(DE_2_down,by="Symbol") %>%
  dplyr::semi_join(DE_3_down,by="Symbol") %>%
  t() %>%
  as.character()-> genes_down_in_all_pairs

c(genes_up_in_all_pairs,genes_down_in_all_pairs) -> genes_DE_in_all_pairs

# noiseq result overlap ---------------------------------------------------
# BGI result filter
# filter condition: FC1.5 probability0.8
DE_all %>%
  dplyr::select(Symbol,Probability,`log2FoldChange(G1/Control)`) %>%
  dplyr::filter(Probability>=0.8) %>%
  dplyr::filter(abs(`log2FoldChange(G1/Control)`)>=fc_threshold) %>%
  dplyr::rename("log2FC"=`log2FoldChange(G1/Control)`) %>%
  dplyr::mutate(`G1/Control`=ifelse(log2FC>0,"Up","Down")) -> DE_all_0.8_1.5

# overlap -----------------------------------------------------------------

DE_all_0.8_1.5 %>%
  dplyr::filter(Symbol %in% genes_DE_in_all_pairs) -> DE_all_0.8_1.5_overlap_allPair



# statistic ---------------------------------------------------------------

DE_all_0.8_1.5_overlap_allPair %>%
  dplyr::select(`G1/Control`)  %>%
  table()

DE_all_0.8_1.5_overlap_allPair %>%
  dplyr::select(Probability,log2FC) %>%
  dplyr::mutate(Probability=log10(1-Probability)) %>%
  

# plot --------------------------------------------------------------------

DE_all %>%
  dplyr::rename("log2FC"=`log2FoldChange(G1/Control)`) %>%
  dplyr::mutate(`G1/Control`=ifelse(log2FC>=(0.585) ,"Up","None")) %>%
  dplyr::mutate(`G1/Control`=ifelse(log2FC<=(-0.585),"Down",`G1/Control`)) %>%
  dplyr::mutate(`G1/Control`=ifelse(Probability>=0.8,`G1/Control`,"None")) %>%
  dplyr::mutate(`G1/Control`=ifelse(is.na(Probability),"None",`G1/Control`)) %>%
  dplyr::select(Symbol,log2FC,Probability,`G1/Control`) %>%
  dplyr::mutate(point=ifelse(Symbol %in% genes_DE_in_all_pairs & `G1/Control` !="None","1","2")) %>%
  dplyr::mutate(Probability=-log10(1-Probability)) %>%
  dplyr::mutate(color=ifelse(`G1/Control`=="Up","red","grey")) %>%
  dplyr::mutate(color=ifelse(`G1/Control`=="Down","blue",color)) -> point_ready

library(ggplot2)
point_ready %>%
  ggplot() +
  geom_point(aes(x=log2FC,y=Probability,color=`G1/Control`)) + #,colour=point_ready$color
  ylab("-Log10(1-Probability)") -> p

ggsave(file.path(out_path,"plot/G1-control.MA.plot.pdf"),p,device = "pdf",width = 3,height = 3)
  
