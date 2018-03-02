library(magrittr)

# path --------------------------------------------------------------------
enrich_path <- c("/project/huff/huff/TKI/result/IM_DS/mRNA_DE/DAVID-enrich")

# load data ---------------------------------------------------------------

GO_bp <- readr::read_tsv(file.path(enrich_path,"GO-BP.txt"))
GO_cc <- readr::read_tsv(file.path(enrich_path,"GO-CC.txt"))
GO_mf <- readr::read_tsv(file.path(enrich_path,"GO-MF.txt"))

# merge data --------------------------------------------------------------

rbind(GO_bp,GO_cc) %>% rbind(GO_mf) -> DAVID_enrich

# plot --------------------------------------------------------------------
library(ggplot2)

DAVID_enrich %>%
  dplyr::arrange(Category,Count,desc(PValue)) %>%
  dplyr::select(Term) -> term_rank

DAVID_enrich %>%
  ggplot(mapping = aes(x=Term,y=Count)) +
  geom_bar(aes(fill=PValue),stat = 'identity') +
  scale_x_discrete(limit=term_rank$Term) +
  coord_flip() +
  facet_grid(~Category,scales = "free_x") +
  theme(
    axis.text = element_text(size=15),
    legend.text = element_text(size=10)
  )-> p;p


GO_bp %>%
  # dplyr::filter(PValue<=0.05) %>%
  dplyr::arrange(Category,Count,desc(PValue)) %>%
  dplyr::select(Term) -> GO_bp.term_rank
GO_bp %>%
  # dplyr::filter(PValue<=0.05) %>%
  ggplot(mapping = aes(x=Term,y=Count)) +
  geom_bar(aes(fill=PValue),stat = 'identity') +
  scale_x_discrete(limit=GO_bp.term_rank$Term) +
  labs(title="GO BP enrichment.") +
  coord_flip() -> GO_bp.p

GO_cc %>%
  dplyr::filter(PValue<=0.05) %>%
  dplyr::arrange(Category,Count,desc(PValue)) %>%
  dplyr::select(Term) -> GO_cc.term_rank
GO_cc %>%
  dplyr::filter(PValue<=0.05) %>%
  ggplot(mapping = aes(x=Term,y=Count)) +
  geom_bar(aes(fill=PValue),stat = 'identity') +
  scale_x_discrete(limit=GO_cc.term_rank$Term) +
  labs(title="GO CC enrichment.") +
  coord_flip() -> GO_cc.p

GO_mf %>%
  dplyr::filter(PValue<=0.05) %>%
  dplyr::arrange(Category,Count,desc(PValue)) %>%
  dplyr::select(Term) -> GO_mf.term_rank
GO_mf %>%
  dplyr::filter(PValue<=0.05) %>%
  ggplot(mapping = aes(x=Term,y=Count)) +
  geom_bar(aes(fill=PValue),stat = 'identity') +
  scale_x_discrete(limit=GO_mf.term_rank$Term) +
  labs(title="GO MF enrichment.") +
  coord_flip() -> GO_mf.p