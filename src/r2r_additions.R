
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
