
# figure out change in count of supported T-I with new associations datasets added in

merge2 = read_tsv('data/merge2.tsv.gz', col_types=cols())

merge2 %>%
  filter(comb_norm >= 0.8 & !is.na(comb_norm)) %>%
  filter(assoc_source != 'OTG' | l2g_share >= 0.5) %>%
  group_by(ti_uid) %>%
  slice(1) %>%
  group_by(assoc_source) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup()

merge2 %>%
  filter(comb_norm >= 0.8 & !is.na(comb_norm)) %>%
  filter(assoc_source == 'Genebass') %>%
  View()

merge2 %>%
  filter(comb_norm >= 0.8 & !is.na(comb_norm)) %>%
  filter(assoc_source == 'OMIM') %>%
  View()

merge2 %>%
  filter(!is.na(comb_norm) & comb_norm >= 0.8) %>%
  filter(assoc_source != 'intOGen' & (l2g_share >= 0.5 | assoc_source != 'OTG')) %>%
  group_by(ti_uid) %>%
  summarize(.groups='keep', max_hcatnum = max(hcatnum)) %>%
  ungroup() %>%
  group_by(max_hcatnum) %>%
  summarize(.groups='keep', n_ti = length(unique(ti_uid)))


merge2 %>%
  filter(!is.na(comb_norm) & comb_norm >= 0.8) %>%
  filter(assoc_source != 'intOGen' & (l2g_share >= 0.5 | assoc_source != 'OTG')) %>%
  group_by(ti_uid) %>%
  summarize(.groups='keep', max_ccatnum = max(ccatnum)) %>%
  ungroup() %>%
  group_by(max_ccatnum) %>%
  summarize(.groups='keep', n_ti = length(unique(ti_uid)))


# old answer:
# merge2 %>%
#   filter(!is.na(comb_norm) & comb_norm >= 0.8) %>%
#   filter(assoc_source != 'intOGen' & (l2g_share >= 0.5 | assoc_source != 'OTG')) %>%
#   group_by(ti_uid) %>%
#   summarize(.groups='keep', max_hcatnum = max(hcatnum)) %>%
#   ungroup() %>%
#   group_by(max_hcatnum) %>%
#   summarize(.groups='keep', n_ti = length(unique(ti_uid)))
# 
# max_hcatnum  n_ti
# <dbl> <int>
#   1           1   474
# 2           2   175
# 3           3   322
# 4           4   103
# 5           5   188






# 
# 
# 
# 
# 
# combined_ti_gwas_sans_omim    = pipeline_best(merge2, phase='combined', basis='ti', associations=c('PICCOLO','OTG','Genebass'), lacking=c('OMIM'), verbose=F)
# combined_ti_omim_sans_gwas    = pipeline_best(merge2, phase='combined', basis='ti', associations=c('OMIM'), lacking=c('PICCOLO','OTG','Genebass'), verbose=F)
# combined_ti_gwas_andalso_omim = pipeline_best(merge2, phase='combined', basis='ti', associations=c('PICCOLO','OTG','Genebass'), andalso=c('OMIM'), verbose=F)
# combined_ti_omim_andalso_gwas = pipeline_best(merge2, phase='combined', basis='ti', associations=c('OMIM'), andalso=c('PICCOLO','OTG','Genebass'), verbose=F)
# 
# 
# advancement_rr(combined_ti_gwas_sans_omim)
# advancement_rr(combined_ti_omim_sans_gwas)
# advancement_rr(combined_ti_omim)
# advancement_rr(combined_ti_gwas_andalso_omim) # 3.93
# advancement_rr(combined_ti_omim_andalso_gwas) # 3.93
# # above RS numbers now agree; check that there are no T-I on which they disagree:
# omim_plus_gwas_ti = unique(combined_ti_omim_andalso_gwas$ti_uid[combined_ti_omim_andalso_gwas$target_status=='genetically supported target'])
# gwas_plus_omim_ti = unique(combined_ti_gwas_andalso_omim$ti_uid[combined_ti_gwas_andalso_omim$target_status=='genetically supported target'])
# setdiff(omim_plus_gwas_ti, gwas_plus_omim_ti)
# setdiff(gwas_plus_omim_ti, omim_plus_gwas_ti)
# 
# # experiment with historical being T-I no longer in active
# pp$hcat[!is.na(pp$acat)] = NA
# pp$hcatnum[!is.na(pp$acat)] = NA
# 
# merge2$hcat[!is.na(merge2$acat)] = NA
# merge2$hcatnum[!is.na(merge2$acat)] = NA
# 
# 
# 
# combined_ti_genebass_misskat  = pipeline_best(merge2, phase='combined', basis='ti', associations='Genebass', genebass_subcat='missenseLC skat', verbose=F)
# combined_ti_genebass_misburd  = pipeline_best(merge2, phase='combined', basis='ti', associations='Genebass', genebass_subcat='missenseLC burden', verbose=F)
# combined_ti_genebass_lofskat  = pipeline_best(merge2, phase='combined', basis='ti', associations='Genebass', genebass_subcat='pLoF skat', verbose=F)
# combined_ti_genebass_lofburd  = pipeline_best(merge2, phase='combined', basis='ti', associations='Genebass', genebass_subcat='pLoF burden', verbose=F)
# 
# advancement_rr(combined_ti_genebass_misskat)
# advancement_rr(combined_ti_genebass_misburd)
# advancement_rr(combined_ti_genebass_lofskat)
# advancement_rr(combined_ti_genebass_lofburd)

