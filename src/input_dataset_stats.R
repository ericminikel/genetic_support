options(stringsAsFactors=F)
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
if(interactive()) {
  setwd('~/d/sci/src/genetic_support')
}

omim_relational_all = read_tsv('../digap/output/omim_relational_all.tsv', col_types=cols())
mendelian_curation = read_tsv('../digap/data/curated/mendelian_curation_corrected.tsv', guess_max = 10000)

omim_relational_all %>%
  left_join(mendelian_curation, by=c('gene_mim','phenotype_mim')) %>%
  filter(real) %>%
  filter(is.na(final_call) | final_call != 'Association not established') %>%
  filter(!duplicated(paste0(gene_mim, phenotype_mim))) -> omim_filtered

nrow(omim_relational_all)
nrow(omim_filtered)


piccolo = read_tsv('../digap/data/piccolo/media-6.txt', col_types=cols()) %>% clean_names()
piccolo$gene_trait = paste0(piccolo$hgnc_idx, '-', piccolo$trait)
length(unique(piccolo$gene_trait))
length(unique(piccolo$gene_trait[piccolo$h4 >= 0.9 & 
                                   piccolo$pval_idx < 1e-5  & 
                                   piccolo$gwas_pvals < 5e-8]))



otg_l2g_v2g = read_tsv('../digap/output/otg2209/l2g.tsv.gz', col_types=cols())
v2d_uniq_assoc = read_tsv('../digap/output/otg2209/v2d_uniq_assoc.tsv.gz', col_types=cols())
uniq_traits = read_tsv('../digap/output/otg2209/uniq_traits.tsv', col_types=cols())
otg_traits = read_tsv('../digap/output/otg2209/traits.tsv.gz', col_types=cols())

nrow(v2d_uniq_assoc)
nrow(uniq_traits)
nrow(otg_traits)

sum(otg_l2g_v2g$l2g_share >= 0.5)

# note that OTG "unique rows" or "unique hits" represent different rows in the original L2G data, which
# are unique on gene - original_trait - lead SNP, though we didn't capture lead SNP in the assoc
# table so instead they are _almost_ unique on gene - original_trait - pval
# and even closer to unique on gene - original_trait - pval - l2g_share




genebass = read_tsv('../digap/hail/genebass_plof_skat_1e-4.tsv', col_types=cols()) %>% 
  clean_names() %>%
  filter(pvalue < 1e-5)
nrow(genebass)



intogen_smry = read_tsv('../digap/output/intogen_smry.tsv', col_types=cols())
nrow(intogen_smry)
