




### version of some panels without genetic insight filter

cat(file=stderr(), 'Generating alternate pipeline tables without genetic insight filter...')

combined_ti_germline_unfiltered = pipeline_best(merge2, phase='combined', basis='ti', require_insight=F, include_missing=F, verbose=F)
combined_ti_omim_unfiltered     = pipeline_best(merge2, phase='combined', basis='ti', require_insight=F, associations=c('OMIM'), verbose=F)
combined_ti_genebass_unfiltered = pipeline_best(merge2, phase='combined', basis='ti', require_insight=F, associations=c('Genebass'), verbose=F)
combined_ti_otg_unfiltered      = pipeline_best(merge2, phase='combined', basis='ti', require_insight=F, associations=c('OTG'), verbose=F)
combined_ti_pic_unfiltered      = pipeline_best(merge2, phase='combined', basis='ti', require_insight=F, associations=c('PICCOLO'), verbose=F)
combined_ti_gwas_unfiltered     = pipeline_best(merge2, phase='combined', basis='ti', require_insight=F, associations=c('PICCOLO','OTG','Genebass'), verbose=F)
combined_ti_gwascat_unfiltered  = pipeline_best(merge2, phase='combined', basis='ti', require_insight=F, otg_subcat='GWAS Catalog', verbose=F)
combined_ti_ukbb_unfiltered     = pipeline_best(merge2, phase='combined', basis='ti', require_insight=F, otg_subcat='Neale UKBB', verbose=F)
combined_ti_finngen_unfiltered  = pipeline_best(merge2, phase='combined', basis='ti', require_insight=F, otg_subcat='FinnGen', verbose=F)
cat(file=stderr(), 'done!\n')

assoc_source_rr_forest_unfiltered = tibble(label=c('All germline','OMIM','All GWAS','All OTG','GWAS Catalog','Neale UKBB','FinnGen','PICCOLO','Genebass'),
                                pipeline_obj=c('combined_ti_germline_unfiltered','combined_ti_omim_unfiltered','combined_ti_gwas_unfiltered','combined_ti_otg_unfiltered','combined_ti_gwascat_unfiltered','combined_ti_ukbb_unfiltered','combined_ti_finngen_unfiltered','combined_ti_pic_unfiltered','combined_ti_genebass_unfiltered')) %>%
  mutate(y=max(row_number()) - row_number() + 1)
for (i in 1:nrow(assoc_source_rr_forest_unfiltered)) {
  pipeline_obj = get(assoc_source_rr_forest_unfiltered$pipeline_obj[i])
  rr_obj = advancement_rr(pipeline_obj)
  assoc_source_rr_forest_unfiltered[i,c('mean','l95','u95')] = rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  assoc_source_rr_forest_unfiltered[i,c('numerator')]        = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  assoc_source_rr_forest_unfiltered[i,c('denominator')]      = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
}

plot_forest(assoc_source_rr_forest_unfiltered, xlims=c(0,5), xstyle='ratio', mar=c(2.5, 8, 3, 8))
mtext(side=1, line=1.6, text='RS')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


assoc_source_rr_forest_unfiltered %>%
  select(label, rs_unfiltered = mean) %>%
  inner_join(select(assoc_source_rr_forest, label, rs_filtered=mean), by = 'label') %>%
  mutate(difference = rs_unfiltered - rs_filtered) %>%
  mutate(y_offset = case_when(label %in% c('All germline', 'All GWAS') ~ -0.03,
                              TRUE ~ 0 )) -> assoc_source_filt_un
par(mar=c(3,3,3,1))
xlims = c(2,4)
ylims = c(2,4)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:4, labels=NA)
axis(side=1, at=0:4, line=-0.5, lwd=0)
mtext(side=1, line=1.6, text='RS filtered for genetic insight')
axis(side=2, at=0:4, labels=NA)
axis(side=2, at=0:4, line=-0.5, lwd=0, las=2)
mtext(side=2, line=1.6, text='RS unfiltered')
abline(a=0, b=1, col='black', lty=3)
points(assoc_source_filt_un$rs_filtered, assoc_source_filt_un$rs_unfiltered, pch=20)
text(assoc_source_filt_un$rs_filtered, assoc_source_filt_un$rs_unfiltered + assoc_source_filt_un$y_offset, labels=assoc_source_filt_un$label, pos=4, cex=0.7)

assoc_source_filt_un %>%
  select(-y_offset) -> assoc_source_filt_un_out
write_supp_table(assoc_source_filt_un_out, 'RS by association source with and without genetic insight filter.')


areas_all %>%
  select(area, topl, filter,color, rs_filtered = rs_mean) -> areas_all_filt_un

for (i in 1:nrow(areas_all)) {
  combined_ti_unfilt_area = subset_by_area(combined_ti_unfiltered, areas_all_filt_un$topl[i], areas_all_filt_un$filter[i])
  area_rr_obj = advancement_rr(combined_ti_unfilt_area)
  areas_all_filt_un[i,c('rs_unfiltered','rs_l95','rs_u95')]         = area_rr_obj[area_rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
}
areas_all_filt_un %>%
  mutate(difference = rs_unfiltered - rs_filtered) %>%
  select(area, rs_filtered, rs_unfiltered, difference) -> areas_all_filt_un_out
write_supp_table(areas_all_filt_un_out, 'RS by therapy area with and without genetic insight filter.')

par(mar=c(3,3,3,1))
xlims = c(0,4)
ylims = c(0,4)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:4, labels=NA)
axis(side=1, at=0:4, line=-0.5, lwd=0)
mtext(side=1, line=1.6, text='RS filtered for genetic insight')
axis(side=2, at=0:4, labels=NA)
axis(side=2, at=0:4, line=-0.5, lwd=0, las=2)
mtext(side=2, line=1.6, text='RS unfiltered')
abline(a=0, b=1, col='black', lty=3)
points(areas_all_filt_un$rs_filtered, areas_all_filt_un$rs_unfiltered, col=areas_all_filt_un$color, pch=20)
text(areas_all_filt_un$rs_filtered, areas_all_filt_un$rs_unfiltered, labels=areas_all_filt_un$area, col=areas_all_filt_un$color, pos=4, cex=0.7)









# figure out change in count of supported T-I with new associations datasets added in

old = read_tsv('ignore/merge2_2023-05.tsv.gz', col_types=cols())
new = read_tsv('data/merge2.tsv.gz', col_types=cols())

merge2 = new

merge2 %>%
  filter(comb_norm >= 0.8 & !is.na(comb_norm)) %>%
  filter(assoc_source != 'OTG' | l2g_share >= 0.5) %>%
  group_by(ti_uid) %>%
  slice(1) %>%
  group_by(assoc_source) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() -> gs_ti_by_source
gs_ti_by_source
sum(gs_ti_by_source$n) # up from 1763 to 1773

merge2 %>%
  filter(comb_norm >= 0.8 & !is.na(comb_norm)) %>%
  filter(assoc_source == 'Genebass') -> genebass_ti
nrow(genebass_ti) # up from 85 to 257

merge2 %>%
  filter(comb_norm >= 0.8 & !is.na(comb_norm)) %>%
  filter(assoc_source == 'OMIM') -> omim_ti
nrow(omim_ti) # down from 484 to 457

new %>%
  filter(comb_norm >= 0.8 & !is.na(comb_norm)) %>%
  filter(assoc_source == 'OMIM') -> omim_ti_new
old %>%
  filter(comb_norm >= 0.8 & !is.na(comb_norm)) %>%
  filter(assoc_source == 'OMIM') -> omim_ti_old

# pairs that lost support in the new version
omim_support_lost = setdiff(omim_ti_old$ti_uid, omim_ti_new$ti_uid)
new %>% filter(ti_uid %in% omim_support_lost & assoc_source=='OMIM') %>% View()
old %>% filter(ti_uid %in% omim_support_lost & assoc_source=='OMIM') %>% View()
# ok it is just slight changes in similarity. e.g. ABCA1-D006949 (hyperlipidemias) was previously supported by D052456 hypoalphalipoproteinemias
# at a similarity of 0.804 but now it has dropped to 0.797
old %>% 
  filter(ti_uid %in% omim_support_lost & assoc_source=='OMIM') %>%
  filter(comb_norm >= 0.8) %>%
  select(ti_uid, assoc_mesh_id, assoc_mesh_term, comb_norm) -> old_matches
new %>% 
  filter(ti_uid %in% omim_support_lost & assoc_source=='OMIM') %>%
  select(ti_uid, assoc_mesh_id, assoc_mesh_term, comb_norm) -> new_matches
old_matches %>%
  left_join(new_matches, by=c('ti_uid','assoc_mesh_id','assoc_mesh_term'), suffix=c('_old','_new')) %>% View()
# ok most are from just above 0.80 to just below.
# however a handful have indication_mesh_id that is an SCR, sometimes with comb_norm==1.0 and the match is now entirely missing in OMIM
# e.g.:
new %>%
  filter(ti_uid=="MMUT-C537358") %>% View()
# this was an _exact_ match on the SCR with sim 1.0, but now OMIM is mapped to D000592 Amino Acid Metabolism, Inborn Errors, which
# although an excellent mapping, has comb_norm 0 because the matrix does poorly at handling SCRs.
# solution: on Pharmaprojects as on OMIM, roll up each SCR to its "mapped to" main heading.

new %>%
  filter(!is.na(comb_norm) & comb_norm >= 0.8) %>%
  filter(assoc_source != 'intOGen' & (l2g_share >= 0.5 | assoc_source != 'OTG')) %>%
  group_by(ti_uid) %>%
  summarize(.groups='keep', max_hcatnum = max(hcatnum)) %>%
  ungroup() %>%
  group_by(max_hcatnum) %>%
  summarize(.groups='keep', n_ti = length(unique(ti_uid)))


new %>%
  filter(!is.na(comb_norm) & comb_norm >= 0.8) %>%
  filter(assoc_source != 'intOGen' & (l2g_share >= 0.5 | assoc_source != 'OTG')) %>%
  group_by(ti_uid) %>%
  summarize(.groups='keep', max_ccatnum = max(ccatnum)) %>%
  ungroup() %>%
  group_by(max_ccatnum) %>%
  summarize(.groups='keep', n_ti = length(unique(ti_uid))) %>%
  ungroup() -> new_gs_ti_by_phase

old %>%
  filter(!is.na(comb_norm) & comb_norm >= 0.8) %>%
  filter(assoc_source != 'intOGen' & (l2g_share >= 0.5 | assoc_source != 'OTG')) %>%
  group_by(ti_uid) %>%
  summarize(.groups='keep', max_ccatnum = max(ccatnum)) %>%
  ungroup() %>%
  group_by(max_ccatnum) %>%
  summarize(.groups='keep', n_ti = length(unique(ti_uid))) %>%
  ungroup() -> old_gs_ti_by_phase

old_gs_ti_by_phase %>%
  inner_join(new_gs_ti_by_phase, by='max_ccatnum', suffix=c('_old','_new')) -> old_v_new_gs_ti_by_phase

sum(old_v_new_gs_ti_by_phase$n_ti_old[old_v_new_gs_ti_by_phase$max_ccatnum > 1])
sum(old_v_new_gs_ti_by_phase$n_ti_new[old_v_new_gs_ti_by_phase$max_ccatnum > 1])

# down from 779 to 743



new %>%
  filter(!is.na(comb_norm) & comb_norm >= 0.8) %>%
  group_by(ti_uid) %>%
  summarize(.groups='keep', max_ccatnum = max(ccatnum)) %>%
  ungroup() %>%
  group_by(max_ccatnum) %>%
  summarize(.groups='keep', n_ti = length(unique(ti_uid))) %>%
  ungroup() -> new_gs_ti_by_phase_unfilt

old %>%
  filter(!is.na(comb_norm) & comb_norm >= 0.8) %>%
  group_by(ti_uid) %>%
  summarize(.groups='keep', max_ccatnum = max(ccatnum)) %>%
  ungroup() %>%
  group_by(max_ccatnum) %>%
  summarize(.groups='keep', n_ti = length(unique(ti_uid))) %>%
  ungroup() -> old_gs_ti_by_phase_unfilt

old_gs_ti_by_phase_unfilt %>%
  inner_join(new_gs_ti_by_phase_unfilt, by='max_ccatnum', suffix=c('_old','_new')) -> old_v_new_gs_ti_by_phase_unfilt

# both are exactly 2153:
sum(old_v_new_gs_ti_by_phase_unfilt$n_ti_old)
sum(old_v_new_gs_ti_by_phase_unfilt$n_ti_new)
# it is a coincidence, because the numbers in each phase actually differ:
old_v_new_gs_ti_by_phase_unfilt
