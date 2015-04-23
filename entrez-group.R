library(dplyr)
library(ggplot2)

# WARNING: do not attempt without a minimum of 32 GB of RAM

# read ucsc gene data
erc.df <- file.path('data', 'erc_mam33.tsv.gz') %>%
  read.delim(colClasses=c('character', 'integer', 'character', 'integer', 'numeric'), na.strings='')

# collapse into entrez genes
entrez.df <- erc.df %>%
  na.omit() %>%
  dplyr::group_by(source_entrez, target_entrez) %>%
  dplyr::summarize(
    correlation = mean(correlation),
    n_ucsc_genes = n()
  ) %>%
  dplyr::ungroup()

# write entrez gene values to disk
write.gz.tsv <- function(df, path) {
  gz <- gzfile(path, 'w')
  write.table(df, gz, sep='\t', quote=FALSE, row.names=FALSE)
  close(gz)
}

write.gz.tsv(entrez.df, file.path('data', 'erc_mam33-entrez.tsv.gz'))

# write gene pairs with ERC >= 0.6
entrez.df %>%
  dplyr::filter(correlation >= 0.6) %>%
  write.gz.tsv(file.path('data', 'erc_mam33-entrez-gt-0.6.tsv.gz'))

# plot ERC values
gg.erc <- entrez.df %>%
  ggplot(aes(x = correlation)) +
  geom_histogram(binwidth=0.01, alpha=0.7) +
  xlab('Evolutionary Rate Covariation Value') + ylab('Count') + 
  theme_bw() + theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

file.path('figure', 'erc-value-dist.png') %>%
  ggsave(gg.erc, width=6, height=3, dpi=200)

# Signed ERC values
sign.df <- entrez.df %>%
  dplyr::mutate(sign = ifelse(correlation < 0, 'negative', 'positive')) %>%
  dplyr::mutate(correlation = abs(correlation))

gg.sign <- sign.df %>%
  ggplot(aes(x = correlation)) +
  geom_histogram(aes(fill = sign), binwidth=0.01, alpha=0.9, position="fill") +
  xlab('Absolute Evolutionary Rate Covariation Value') + ylab('Probability') + 
  theme_bw() + theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) + 
  theme(legend.justification=c(1, 1), legend.position=c(1, 1))

file.path('figure', 'erc-signed-dist.png') %>%
  ggsave(gg.sign, width=6, height=3, dpi=200)

# plot UCSC genes per Entrez gene
gg.map <- entrez.df %>%
  ggplot(aes(x = n_ucsc_genes)) +
  geom_histogram(binwidth=1, alpha=0.7) +
  xlab('UCSC Genes per Entrez Gene') + ylab('Count') + 
  theme_bw() + theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

file.path('figure', 'ucsc-per-entrez-dist.png') %>%
  ggsave(gg.map, width=6, height=3, dpi=200)


# entrez to ucsc mapping
map.df <- dplyr::bind_rows(
  dplyr::transmute(erc.df, ucsc_gene = source_ucsc, entrez_gene = source_entrez),
  dplyr::transmute(erc.df, ucsc_gene = target_ucsc, entrez_gene = target_entrez)
  ) %>%
  dplyr::filter(!is.na(entrez_gene)) %>%
  dplyr::distinct()

map.df %>%
  write.table('data/ucsc-to-entrez-gene.tsv', sep='\t', quote=FALSE, row.names=FALSE)

