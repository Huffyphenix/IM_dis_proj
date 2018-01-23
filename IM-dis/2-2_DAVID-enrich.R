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

rbind(GO_cc,GO_mf) %>% rbind(kegg) -> DAVID_enrich


# plot --------------------------------------------------------------------
GO_bp %>%
  dplyr::filter(PValue<=0.05) %>%
  dplyr::arrange(Category,Count,desc(PValue)) %>%
  dplyr::select(Term) -> GO_bp.term_rank
GO_bp %>%
  dplyr::filter(PValue<=0.05) %>%
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

kegg %>%
  dplyr::filter(PValue<=0.05) %>%
  dplyr::arrange(Category,Count,desc(PValue)) %>%
  dplyr::select(Term) -> kegg.term_rank
kegg %>%
  dplyr::filter(PValue<=0.05) %>%
  ggplot(mapping = aes(x=Term,y=Count)) +
  geom_bar(aes(fill=PValue),stat = 'identity') +
  scale_x_discrete(limit=kegg.term_rank$Term) +
  labs(title="KEGG pathway enrichment.") +
  coord_flip() -> kegg.p



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


DAVID_enrich %>%
  dplyr::arrange(Category,Count,desc(PValue)) %>%
  dplyr::select(Term) -> term_rank

DAVID_enrich %>%
  ggplot(mapping = aes(x=Term,y=Count)) +
  geom_bar(aes(fill=PValue),stat = 'identity') +
  scale_x_discrete(limit=term_rank$Term) +
  coord_flip() +
  facet_grid(~Category,scales = "free_x") -> p


multiplot(GO_bp.p,p,cols=2)
ggsave(filename = file.path(out_path,"DAVID_enrich.pdf"),device = "pdf",width = 12,height = 6)
