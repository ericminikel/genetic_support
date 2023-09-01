start_time = Sys.time()
cat(file=stderr(), 'Loading dependencies...')

options(stringsAsFactors=F)
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))
if(interactive()) {
  setwd('~/d/sci/src/genetic_support')
}

cat(file=stderr(), '\nReading in data...')

best_names = read_tsv('../digap/output/mesh_best_names.tsv.gz', col_types=cols())
mesh_tree = read_delim('../digap_pipeline/data/mesh_tree.tsv', col_names = c("id", "tree"), col_types=cols())
remaps = read_tsv('data_received/mesh_remappings.tsv', col_types=cols())

meta_hcat = read_tsv('../digap_pipeline/data/meta_hcat.tsv', col_types=cols())
meta_acat = read_tsv('../digap_pipeline/data/meta_acat.tsv', col_types=cols())
meta_ccat = read_tsv('../digap_pipeline/data/meta_ccat.tsv', col_types=cols())

areas = read_tsv('data/areas.tsv', col_types=cols())
universe = read_tsv('../digap/data/universe/universe.tsv', col_types=cols())

assoc = read_tsv('../digap/output/all_mesh_assoc.tsv', col_types=cols())

pp_new = read.xlsx('data_received/target_indication_pair_with_all_info.xlsx', na.strings = 'NA') %>%
  clean_names() %>%
  mutate(highest_status_reached = na_if(highest_status_reached,'NA')) %>%
  mutate(current_status = na_if(current_status,'NA')) %>%
  mutate(highest_status_reached_single_target = na_if(highest_status_reached_single_target,'NA')) %>%
  mutate(current_status_single_target = na_if(current_status_single_target,'NA')) %>%
  select(gene,
         indication_mesh_id = mesh_id,
         highest_status_reached,
         current_status,
         highest_status_reached_single_target,
         current_status_single_target,
         orphan=is_orphan_drug,
         year_launch) %>%
  left_join(remaps, by=c('indication_mesh_id'='old')) %>%
  mutate(indication_mesh_id = case_when(indication_mesh_id %in% remaps$old ~ new,
                                        TRUE ~ indication_mesh_id))

sim = read_tsv('data/sim.tsv.gz', col_types=cols())

####
# Pharmaprojects analytic dataset
####

cat(file=stderr(), 'done!\nPreparing Pharmaprojects...')

minimum_date = as.Date('2000-01-01')

nany = function(x) {
  case_when(any(as.logical(x)) ~ T,
            any(as.logical(x) %in% F) ~ F,
            all(is.na(as.logical(x))) ~ NA)
}


pp_new %>% 
  mutate(ti_uid = paste0(gene,'-',indication_mesh_id)) %>%
  mutate(hcat = factor(highest_status_reached, levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'), ordered=T)) %>%
  mutate(acat = factor(current_status, levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'), ordered=T)) %>%
  mutate(ccat = pmax(hcat,acat,na.rm=T)) %>%
  filter(!is.na(hcat)) %>%
  mutate(oto_hcat = factor(highest_status_reached_single_target, levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'), ordered=T)) %>%
  mutate(oto_acat = factor(current_status_single_target, levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'), ordered=T)) %>%
  mutate(oto_ccat = pmax(oto_hcat,oto_acat,na.rm=T)) %>%
  mutate(succ_p_1 = case_when(ccat %in% c('Phase I','Phase II','Phase III','Launched') ~ TRUE,
                              hcat %in% c('Preclinical') ~ FALSE,
                              acat %in% c('Preclinical') ~ NA,
                              TRUE ~ NA),
         succ_1_2 = case_when(ccat %in% c('Phase II','Phase III','Launched') ~ TRUE,
                              hcat %in% c('Phase I') ~ FALSE,
                              acat %in% c('Phase I') ~ NA,
                              TRUE ~ NA),
         succ_2_3 = case_when(ccat %in% c('Phase III','Launched') ~ TRUE,
                              hcat %in% c('Phase II') ~ FALSE,
                              acat %in% c('Phase II') ~ NA,
                              TRUE ~ NA),
         succ_3_a = case_when(ccat %in% c('Launched') ~ TRUE,
                              hcat %in% c('Phase III') ~ FALSE,
                              acat %in% c('Phase III') ~ NA,
                              TRUE ~ NA)) %>%
  left_join(best_names, by=c('indication_mesh_id'='id')) %>%
  rename(indication_mesh_term = labeltext) %>%
  filter(!grepl('iagnos',indication_mesh_term)) %>% # remove approvals for diagnostic applications (as opposed to treatment)
  left_join(meta_acat, by=c('acat'='cat')) %>%
  rename(acatnum = num) %>%
  left_join(meta_hcat, by=c('hcat'='cat')) %>%
  rename(hcatnum = num) %>%
  left_join(meta_ccat, by=c('ccat'='cat')) %>%
  rename(ccatnum = num) %>%
  select(ti_uid, gene, indication_mesh_id, indication_mesh_term,
         hcat, acat, ccat, hcatnum, acatnum, ccatnum,
         succ_p_1, succ_1_2, succ_2_3, succ_3_a,
         oto_hcat, oto_acat, oto_ccat, orphan, year_launch) -> pp_analytic

write_tsv(pp_analytic, 'data/pp.tsv', na='')

pp_launched_alltime = read.xlsx('nonrelease/input_sources-2023-01-20/target_indication_highest_status_all_years.xlsx') %>%
  rename(indication_mesh_id = mesh_id) %>%
  mutate(ti_uid = paste0(gene,'-',indication_mesh_id)) %>%
  mutate(hcat = factor(highest_status_reached, levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'), ordered=T)) %>%
  left_join(best_names, by=c('indication_mesh_id'='id')) %>%
  rename(indication_mesh_term = labeltext) %>%
  filter(!grepl('iagnos',indication_mesh_term)) %>% 
  filter(!is.na(hcat) & !is.na(gene) & !is.na(indication_mesh_id)) %>%
  filter(hcat=='Launched') %>%
  select(ti_uid, gene, indication_mesh_id, indication_mesh_term) 

write_tsv(pp_launched_alltime, 'data/pp_launched_alltime.tsv', na='')

pp_drugs = read.xlsx('data_received/drug_target_mesh_after_2000_202212222.xlsx') %>%
  clean_names() %>% as_tibble()

meta_acat %>%
  inner_join(pp_drugs, by=c('cat'='highest_status_reached')) %>%
  group_by(drug_id) %>%
  summarize(.groups='keep', maxphasenum=max(num)) %>%
  ungroup() %>%
  inner_join(meta_acat, by=c('maxphasenum'='num')) %>%
  rename(maxphase = cat) -> drugmaxphase

pp_drugs %>%
  inner_join(meta_acat, by=c('highest_status_reached'='cat')) %>%
  rename(phasenum = num,
         phase=highest_status_reached) %>%
  inner_join(drugmaxphase, by=c('drug_id')) %>%
  mutate(ti_uid = paste0(gene, '-',mesh_id)) %>%
  group_by(ti_uid, phasenum, phase, maxphasenum, maxphase) %>%
  summarize(.groups='keep',
            n_drugs = length(unique(drug_id))) -> drug_phase_summary

write_tsv(drug_phase_summary, 'data/drug_phase_summary.tsv', na='')

pp_drugs %>%
  group_by(mesh_id) %>%
  summarize(.groups='keep', n_drugs=length(unique(drug_id))) -> drugs_per_indic

write_tsv(drugs_per_indic, 'data/drugs_per_indic.tsv', na='')


####
# Genetic associations analytic dataset
####

cat(file=stderr(), 'done!\nPreparing genetic associations...')


assoc %>%
  filter(!is.na(gene)) %>% # no point in keeping gene-less assocs 
  filter(!is.na(mesh_id)) %>% # no point in keeping mesh-unmapped assocs
  filter(!(source=='OTG' & is.na(l2g_share))) %>% # for OTG only accept entries with L2G
  mutate(extra_info = ifelse(source=='OMIM',NA,extra_info)) %>% # remove disease mechanism curation info
  select(-mesh_term) %>%
  left_join(best_names, by=c('mesh_id'='id')) %>% # map latest/best names for each mesh id
  rename(mesh_term = labeltext) %>%
  select(-v2g_share, -v2g_rank) -> assoc_analytic

write_tsv(assoc_analytic, 'data/assoc.tsv.gz', na='')

####
# Pre-staged merged data (Pharmaprojects x genetic associations)
####

cat(file=stderr(), 'done!\nStaging merged data...')

pp_analytic %>%
  inner_join(assoc_analytic, by='gene', suffix=c('','_a')) %>%
  select(ti_uid:succ_3_a, orphan, year_launch,
         arow,
         assoc_mesh_id=mesh_id,
         assoc_mesh_term=mesh_term,
         assoc_source=source,
         assoc_info=extra_info,
         original_trait,
         original_link,
         assoc_year=year,
         pic_qtl_pval=pic_qtl_pval,
         pic_h4=pic_h4,
         af_gnomad_nfe=af_gnomad_nfe,
         l2g_share=l2g_share,
         l2g_rank=l2g_rank) -> merge1_present

pp_analytic %>%
  filter(!is.na(gene)) %>% 
  select(ti_uid:succ_3_a, orphan, year_launch) %>%
  mutate(arow=as.integer(NA),
         assoc_mesh_id=as.character(NA),
         assoc_mesh_term=as.character(NA),
         assoc_source=as.character(NA),
         assoc_info=as.character(NA),
         original_trait=as.character(NA),
         original_link=as.character(NA),
         assoc_year=as.integer(NA),
         pic_qtl_pval=as.numeric(NA),
         pic_h4=as.numeric(NA),
         af_gnomad_nfe=as.numeric(NA),
         l2g_share=as.numeric(NA),
         l2g_rank=as.integer(NA)) -> merge1_blanks

merge1 = rbind(merge1_present, merge1_blanks)

merge1 %>%
  left_join(sim, by=c('indication_mesh_id' = 'meshcode_a', 'assoc_mesh_id' = 'meshcode_b')) -> merge2

# the sim matrix does not contain self rows, so manually assign 1
merge2$comb_norm[merge2$indication_mesh_id==merge2$assoc_mesh_id] = 1
# and, assign 0 where either similarity was missing or there was no association or things unspecified:
merge2$comb_norm[is.na(merge2$comb_norm)] = 0
# these should never arise due to filters above:
# sum(is.na(merge2$gene)) # merge2$comb_norm[is.na(merge2$gene)] = 0
# sum(is.na(merge2$indication_mesh_id)) # merge2$comb_norm[merge2$indication_mesh_id==''] = 0

write_tsv(merge2, 'data/merge2.tsv.gz')

####
# Indications - genetic insight, areas, and match
####

cat(file=stderr(), 'done!\nCharacterizing indications...')

# indications table
v2d_uniq_assoc = read_tsv('../digap/output/otg2209/v2d_uniq_assoc.tsv.gz',col_types=cols())
v2d_uniq_assoc$k1m = formatC(round(v2d_uniq_assoc$lead_pos/1e6)*1e6,format='d',width=9,flag='0')
v2d_uniq_assoc$locus = paste0(v2d_uniq_assoc$lead_chrom,':',v2d_uniq_assoc$k1m)
assoc_analytic %>%
  group_by(original_trait, mesh_id) %>%
  summarize(.groups='keep') %>%
  arrange(original_trait, mesh_id) -> assoc_map
v2d_uniq_assoc$mesh_id = assoc_map$mesh_id[match(v2d_uniq_assoc$trait_reported,assoc_map$original_trait)]

pp_analytic %>%
  group_by(indication_mesh_id, indication_mesh_term) %>%
  summarize(.groups='keep', 
            n_targets = length(unique(gene)),
            maxacatnum = max(acatnum),
            maxccatnum = max(ccatnum)) %>%
  ungroup() %>%
  arrange(desc(maxccatnum), desc(n_targets)) -> pp_indications1

# # alternative analysis including all genetically associated mesh IDs
# assoc %>%
#   filter(!(mesh_id %in% pp_indications1$indication_mesh_id)) %>%
#   select(indication_mesh_id = mesh_id, indication_mesh_term = mesh_term) %>%
#   distinct(indication_mesh_id, indication_mesh_term) %>%
#   mutate(n_targets = 0, maxacatnum = NA, maxccatnum = NA) -> pp_indications2
# rbind(pp_indications1, pp_indications2) -> pp_indications
pp_indications1 -> pp_indications

sim %>%
  filter(comb_norm >= 0.8,
         meshcode_a %in% pp_analytic$indication_mesh_id) -> mesh_sim_to_join

assoc %>%
  filter(source %in% c('OMIM','TCGA')) %>% 
  inner_join(mesh_sim_to_join, by=c('mesh_id' = 'meshcode_b')) %>%
  inner_join(pp_analytic, by=c('meshcode_a' = 'indication_mesh_id'), keep=T, suffix=c('_a','_p')) %>%
  group_by(indication_mesh_id) %>%
  summarize(.groups='keep', n_omim_genes=n_distinct(gene_a)) %>%
  filter(n_omim_genes >= 1) -> omimtcga_insight

v2d_uniq_assoc %>%
  filter(!is.na(locus)) %>%
  group_by(mesh_id, locus) %>%
  summarize(.groups='keep') -> gwas_loci

pp_analytic %>%
  filter(indication_mesh_id %in% mesh_sim_to_join$meshcode_a) -> pp_to_join

gwas_loci %>%
  filter(mesh_id %in% mesh_sim_to_join$meshcode_b) %>%
  inner_join(mesh_sim_to_join, by=c('mesh_id' = 'meshcode_b')) %>%
  inner_join(pp_to_join, by=c('meshcode_a' = 'indication_mesh_id'), keep=T) %>%
  group_by(indication_mesh_id) %>%
  summarize(.groups='keep', n_gwas_loci=n_distinct(locus)) %>%
  filter(n_gwas_loci >= 1) -> gwas_insight

pp_indications$n_omimtcga_genes = omimtcga_insight$n_omim_genes[match(pp_indications$indication_mesh_id, omimtcga_insight$indication_mesh_id)]
pp_indications$n_omimtcga_genes[is.na(pp_indications$n_omimtcga_genes)] = 0
pp_indications$n_gwas_loci = gwas_insight$n_gwas_loci[match(pp_indications$indication_mesh_id, gwas_insight$indication_mesh_id)]
pp_indications$n_gwas_loci[is.na(pp_indications$n_gwas_loci)] = 0
pp_indications$genetic_insight = 'none'
pp_indications$genetic_insight[pp_indications$n_omimtcga_genes >= 1] = 'omim/tcga'
pp_indications$genetic_insight[pp_indications$n_gwas_loci >= 3] = 'gwas'
pp_indications$genetic_insight[pp_indications$n_omimtcga_genes >= 1 & pp_indications$n_gwas_loci >= 3] = 'both'
pp_indications$genetic_insight[is.na(pp_indications$genetic_insight)] = 'none'

# require presence in sim matrix
pp_indications$genetic_insight[!(pp_indications$indication_mesh_id %in% sim$meshcode_a)] = 'none'

mesh_tree$topl = substr(mesh_tree$tree,1,3)
mesh_tree$letter = substr(mesh_tree$tree,1,1)

# first create the general match table without filters
pp_indications %>%
  left_join(mesh_tree, by=c('indication_mesh_id' = 'id')) %>%
  group_by(indication_mesh_id, indication_mesh_term, topl) %>%
  summarize(.groups='keep') %>%
  arrange(indication_mesh_id, topl) -> indic_topl_match_raw
indic_topl_match_raw$topl[is.na(indic_topl_match_raw$topl)] = 'OTH'

topl_maps = read_tsv('data/mesh_2023_topl_maps.tsv',col_types=cols())

indic_topl_match = indic_topl_match_raw
keeps = indic_topl_match$topl %in% topl_maps$topl[topl_maps$topl == topl_maps$map_to]
potential_remaps = indic_topl_match$topl %in% topl_maps$topl[topl_maps$map_to != topl_maps$topl]
# before remapping, if the indication maps to a kept area, simply use that - delete other mappings
to_delete = potential_remaps & indic_topl_match$indication_mesh_id %in% indic_topl_match$indication_mesh_id[keeps]
if (sum(to_delete) > 0) {
  indic_topl_match = indic_topl_match[-which(to_delete),]
}
# now apply remapping
remaps = indic_topl_match$topl %in% topl_maps$topl[topl_maps$map_to != topl_maps$topl]
indic_topl_match$topl[remaps] = topl_maps$map_to[match(indic_topl_match$topl[remaps], topl_maps$topl)]
# bake in the non-cancer filter: delete all rows that are non-C04 but DO have a match to a thing in C04
delete_cancer = indic_topl_match$topl != 'C04' & indic_topl_match$indication_mesh_id %in% indic_topl_match$indication_mesh_id[indic_topl_match$topl == 'C04']
if (sum(delete_cancer) > 0) {
  indic_topl_match = indic_topl_match[-which(delete_cancer),]
}
# bake in the other filter: delete all rows that are OTH where the indication does match a non-OTH area
delete_other = indic_topl_match$topl == 'OTH' & indic_topl_match$indication_mesh_id %in% indic_topl_match$indication_mesh_id[indic_topl_match$topl != 'OTH']
if (sum(delete_other) > 0) {
  indic_topl_match = indic_topl_match[-which(delete_other),]
}
# remove dups (some that mapped to >1 other area have multiple OTH rows)
indic_topl_match %>%
  group_by(indication_mesh_id, indication_mesh_term, topl) %>%
  slice(1) %>%
  ungroup() -> indic_topl_match

# check that the filters worked
# indic_topl_match %>%
#   group_by(indication_mesh_id) %>%
#   summarize(.groups='keep', n_cancer = sum(topl=='C04'), n_non_cancer = sum(topl != 'C04')) %>%
#   ungroup() %>% filter(n_cancer > 0 & n_non_cancer > 0)
# indic_topl_match %>%
#   group_by(indication_mesh_id) %>%
#   summarize(.groups='keep', n_other = sum(topl=='OTH'), n_non_other = sum(topl != 'OTH')) %>%
#   ungroup() %>% filter(n_other > 0 & n_non_other > 0)


indic_topl_match %>%
  inner_join(areas, by='topl') %>%
  group_by(indication_mesh_id) %>%
  summarize(.groups='keep', areas1 = toString(unique(area))) %>%
  arrange(indication_mesh_id) -> indic_areas
pp_indications$areas = indic_areas$areas1[match(pp_indications$indication_mesh_id, indic_areas$indication_mesh_id)]

write_tsv(pp_indications, 'data/indic.tsv', na='')
write_tsv(indic_topl_match, 'data/indic_topl_match.tsv', na='')

####
# T2D analysis
####

pp_dia_raw = read.xlsx('data_received/type_ii_03292023.xlsx')
t2d_gensup = read_tsv('data/t2d/t2d_genetic_support.tsv', col_types=cols())

curated = read_tsv('data/t2d/drug_gene_curation_results.tsv', col_types=cols())
pp_dia_raw$remap = curated$remap[match(pp_dia_raw$gene, curated$gene)]
pp_dia_raw$gene[!is.na(pp_dia_raw$remap)] = pp_dia_raw$remap[!is.na(pp_dia_raw$remap)]


pp_dia_raw %>%
  filter(!is.na(gene)) %>%
  left_join(t2d_gensup, by='gene') %>%
  mutate(gensup_rank = case_when(gensup_omim ~ 1, 
                                 gensup_estab ~ 2,
                                 gensup_novel ~ 3,
                                 TRUE ~ 4)) %>% # where multiple targets, take the one with the "earliest" genetic support
  arrange(gensup_rank) %>%
  group_by(drug_id, drug_primary_name, mesh_id) %>%
  slice(1) %>%
  inner_join(meta_ccat, by=c('highest_status_reached'='cat')) %>%
  rename(indication_mesh_id = mesh_id, indication_mesh_term = mesh_name) %>%
  group_by(gene) %>%
  arrange(desc(num)) %>%
  slice(1) %>%
  ungroup() %>%
  inner_join(meta_ccat, by=c('num')) %>%
  rename(maxphase = cat, maxphasenum=num) %>%
  select(gene, maxphasenum, maxphase, indication_mesh_id, indication_mesh_term) -> pp_dia_processed

write_tsv(pp_dia_processed, 'data/t2d/pp_diabetes.tsv', na='')

####
# Unmodified supporting tables
####

cat(file=stderr(), 'done!\nWriting out additional supporting tables...')

write_tsv(meta_hcat, 'data/meta_hcat.tsv', na='')
write_tsv(meta_acat, 'data/meta_acat.tsv', na='')
write_tsv(meta_ccat, 'data/meta_ccat.tsv', na='')
write_tsv(universe, 'data/universe.tsv', na='')
write_tsv(sim, 'data/sim.tsv.gz', na='')

time_elapsed = Sys.time() - start_time
cat(file=stderr(), 'done!\nAll tasks complete in',round(time_elapsed,1),units(time_elapsed),'.\n')
