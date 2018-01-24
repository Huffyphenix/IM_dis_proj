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

# # G2_1/G1_1 DE genes
# G2_1 %>%
#   dplyr::inner_join(G1_1,by="Symbol") %>%
#   dplyr::mutate(FC = log2(G2_1/G1_1)) %>%
#   dplyr::filter(FC>=fc_threshold) -> DE_1_up
# G2_1 %>%
#   dplyr::inner_join(G1_1,by="Symbol") %>%
#   dplyr::mutate(FC = log2(G2_1/G1_1)) %>%
#   dplyr::filter(FC <= (-fc_threshold)) -> DE_1_down
# 
# # G2_2/G1_2 De genes
# G2_2 %>%
#   dplyr::inner_join(G1_2,by="Symbol") %>%
#   dplyr::mutate(FC = log2(G2_2/G1_2)) %>%
#   dplyr::filter(FC>=fc_threshold) -> DE_2_up
# G2_2 %>%
#   dplyr::inner_join(G1_2,by="Symbol") %>%
#   dplyr::mutate(FC = log2(G2_2/G1_2)) %>%
#   dplyr::filter(FC<=(-fc_threshold)) -> DE_2_down
# 
# # G2_3/G1_3 DE genes
# G2_3 %>%
#   dplyr::inner_join(G1_3,by="Symbol") %>%
#   dplyr::mutate(FC = log2(G2_3/G1_3)) %>%
#   dplyr::filter(FC <= (-fc_threshold)) -> DE_3_down
# G2_3 %>%
#   dplyr::inner_join(G1_3,by="Symbol") %>%
#   dplyr::mutate(FC = log2(G2_3/G1_3)) %>%
#   dplyr::filter(FC >= fc_threshold) -> DE_3_up

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

# at least 2 samples
BGI_DE_all_0.8_1.5 %>%
  dplyr::inner_join(Down_in_2_samples,by="Symbol") -> Down_in_2_and_all_test
BGI_DE_all_0.8_1.5 %>%
  dplyr::inner_join(Up_in_2_samples,by="Symbol") -> Up_in_2_and_all_test
rbind(Down_in_2_and_all_test,Up_in_2_and_all_test) %>%
  readr::write_tsv(path = file.path(out_path,"DE_in_2-AND-in_BGI_0.8_1.5"))

