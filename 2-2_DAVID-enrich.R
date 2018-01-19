library(magrittr)

# path --------------------------------------------------------------------

data_path<- c("/project/huff/huff/TKI/result/mRNA_DE/IM_dis/pathway_enrich")
out_path <- c(file.path(data_path,"plot"))


# load data ---------------------------------------------------------------

GO_bp <- readr::read_tsv(file.path(data_path,"DAVID-GO-BP.txt"))
GO_cc <- readr::read_tsv(file.path(data_path,"DAVID-GO-CC.txt"))
GO_mf <- readr::read_tsv(file.path(data_path,"DAVID-GO-MF.txt"))
kegg <- readr::read_tsv(file.path(data_path,"DAVID-KEGG.txt"))


# merge data --------------------------------------------------------------

rbind(GO_bp,GO_cc) %>% rbind(GO_mf) %>% rbind(kegg) -> DAVID_enrich


# plot --------------------------------------------------------------------
DAVID_enrich %>%
  dplyr::arrange(Category,Count,desc(PValue)) %>%
  dplyr::select(Term) -> term_rank

DAVID_enrich %>%
  ggplot(mapping = aes(x=Term,y=Count)) +
  geom_bar(aes(fill=PValue),stat = 'identity') +
  scale_x_discrete(limit=term_rank$Term) +
  coord_flip() +
  facet_grid(~Category,scales = "free_x") -> p

ggsave(filename = file.path(out_path,"DAVID_enrich.pdf"),device = "pdf",width = 12,height = 6)
