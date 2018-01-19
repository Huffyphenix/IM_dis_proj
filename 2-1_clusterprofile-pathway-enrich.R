library(magrittr)

# path --------------------------------------------------------------------

out_path <- c("/project/huff/huff/TKI/result/mRNA_DE/IM_dis/pathway_enrich")
data_path <- c("/project/huff/huff/TKI/result/mRNA_DE/IM_dis")

# load data ---------------------------------------------------------------

DE_mRNA <- readr::read_tsv(file.path(data_path,"IM_dis_DE_mRNA"))

# get gene list
DE_mRNA %>%
  dplyr::select(Symbol) %>%
  t() %>%
  as.character() -> gene_list
gene.df <- bitr(gene_list, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

# Go enrichment -----------------------------------------------------------

ego_CC <- enrichGO(gene          = gene.df$ENTREZID,
                   # keyType       = "SYMBOL",
                # universe      = names(gene_list),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego_CC)
ego_BP <- enrichGO(gene          = gene.df$ENTREZID,
                   # universe      = names(gene_list),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   # qvalueCutoff  = 0.1,
                   readable      = TRUE)

ego_MF <- enrichGO(gene          = gene.df$ENTREZID,
                   # universe      = names(gene_list),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# gsea --------------------------------------------------------------------

ego3 <- gseGO(geneList     = gene.df$ENTREZID,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)


# DAVID -------------------------------------------------------------------

david <- enrichDAVID(gene = gene.df$ENTREZID,
                     idType = "ENTREZ_GENE_ID",
                     annotation = "KEGG_PATHWAY",
                     david.user = "clusterProfiler@hku.hk")
