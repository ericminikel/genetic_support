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

best_names = read_tsv('data/mesh_best_names.tsv.gz', col_types=cols())
mesh_tree = read_delim('data/mesh_tree.tsv', col_names = c("id", "tree"), col_types=cols())
remaps = read_tsv('data/mesh_remappings.tsv', col_types=cols())
scr_best = read_tsv('data/mesh_scr_to_best.tsv.gz', col_types=cols())

meta_hcat = read_tsv('data/meta_hcat.tsv', col_types=cols())
meta_acat = read_tsv('data/meta_acat.tsv', col_types=cols())
meta_ccat = read_tsv('data/meta_ccat.tsv', col_types=cols())

areas = read_tsv('data/areas.tsv', col_types=cols())
universe = read_tsv('data/universe.tsv', col_types=cols())

assoc = read_tsv('../digap/output/all_mesh_assoc.tsv', col_types=cols())

pp_new = read.xlsx('data_received/target_indication_pair_with_all_info.xlsx', na.strings = 'NA') %>%
  clean_names() %>%
  left_join(scr_best, by=c('mesh_id'='scr')) %>% # join to mapping of SCRs to main headings. note that sum(duplicated(scr_best$scr))==0 so there are no many-to-one joins here
  mutate(mapped_mesh = coalesce(main, mesh_id)) %>% # choose the main that the SCR is mapped to, but if not an SCR, just use the original
  mutate(highest_status_reached = na_if(highest_status_reached,'NA')) %>%
  mutate(current_status = na_if(current_status,'NA')) %>%
  mutate(highest_status_reached_single_target = na_if(highest_status_reached_single_target,'NA')) %>%
  mutate(current_status_single_target = na_if(current_status_single_target,'NA')) %>%
  select(gene,
         indication_mesh_id = mapped_mesh,
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
  as_tibble() %>%
  mutate(ti_uid = paste0(gene,'-',indication_mesh_id)) %>%
  mutate(hcat = factor(ifelse(is.na(current_status) | highest_status_reached=='Launched', highest_status_reached, as.character(NA)), levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'), ordered=T)) %>%
  mutate(acat = factor(ifelse(highest_status_reached=='Launched', 'Launched', current_status), levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'), ordered=T)) %>%
  mutate(ccat = pmax(hcat,acat,na.rm=T)) %>%
  filter(!is.na(ccat)) %>% # only include when either hcat or acat present
  mutate(oto_hcat = factor(ifelse(is.na(current_status_single_target) | highest_status_reached_single_target=='Launched', highest_status_reached_single_target, as.character(NA)), levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'), ordered=T)) %>%
  mutate(oto_acat = factor(ifelse(highest_status_reached_single_target=='Launched', 'Launched', current_status_single_target), levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'), ordered=T)) %>%
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

suzuki_st4 = read_tsv('data/t2d/suzuki_2023_st4.tsv', col_types=cols()) %>% 
  clean_names()
while(any(is.na(suzuki_st4$chromosome))) {
  suzuki_st4 %>%
    mutate(chromosome = case_when(!is.na(chromosome) ~ chromosome,
                                  is.na(chromosome) ~ lag(chromosome))) -> suzuki_st4
}

# st4 %>%
#   mutate(otglink = paste0('https://genetics.opentargets.org/variant/',chromosome,'_',position_bp_b37,'_',risk,'_',other)) -> st4
# write_tsv(st4[,'index_snv'], 'data/t2d/suzuki_lead_rsids.tsv', col_names = F)
# # retrieve b38 matches from OTG variant-index
# -F = fixed, -w = whole word,  -f = from a file
# you'd think you want -m 1 = 1 match only, however, grep interprets this as 1 match _total_ and not _per line of input file_
# instead just do this, which takes several hours:
# zgrep -F -w -f data/t2d/suzuki_lead_rsids.tsv ~/d/sci/src/digap/output/otg_vi_all.tsv.gz > data/t2d/suzuki_lead_b38.tsv
vi_suzuki = read_tsv('data/t2d/suzuki_lead_b38.tsv', col_names=c('pos_id','rsid','af_gnomad_nfe'), col_types=cols())
# write_tsv(vi_suzuki[,'pos_id'], 'data/t2d/suzuki_lead_pos_id.tsv', col_names=F)
# now go back and parse l2g for these. this takes just 1 minute:
# Rscript src/otg_l2g_parser.R data/opentargetsgenetics/2209/l2g/ ~/d/sci/src/genetic_support/data/t2d/suzuki_lead_pos_id.tsv ~/d/sci/src/genetic_support/data/t2d/suzuki_lead_l2g.tsv
l2g_suzuki = read_tsv('data/t2d/suzuki_lead_l2g.tsv', col_types=cols())
l2g_suzuki %>%
  select(pos_id, gene, l2g_share, rank) %>%
  inner_join(vi_suzuki, by='pos_id') %>%
  select(-af_gnomad_nfe) %>%
  inner_join(suzuki_st4, by=c('rsid'='index_snv')) -> st4_with_l2g
# almost all are in variant-index
# > mean(st4$index_snv %in% vi_suzuki$rsid)
# [1] 0.9937936
# but only 34% have an L2G, i.e. only those that are already a hit in a different GWAS.
# note that L2G is assigned per GWAS-SNP-gene combination, but for any SNP-gene it varies very little (1-2%) between
# GWAS, so it should be fine to use the L2G values assigned for other studies for these SNPs
# > mean(vi_suzuki$pos_id %in% l2g_suzuki$pos_id)
# [1] 0.3425712
l2g_suzuki %>%
  filter(l2g_share > 0.5) %>%
  group_by(pos_id, gene) %>%
  summarize(.groups='keep', l2g_share = max(l2g_share)) %>%
  ungroup() %>%
  inner_join(vi_suzuki, by='pos_id') %>%
  select(-af_gnomad_nfe) -> l2g_half
suzuki_st4$l2g_gene = l2g_half$gene[match(suzuki_st4$index_snv, l2g_half$rsid)]
suzuki_st4$pos_id_b38 = vi_suzuki$pos_id[match(suzuki_st4$index_snv, vi_suzuki$rsid)]
#sum(!is.na(suzuki_st4$l2g_gene)) # 346
#mean(!is.na(suzuki_st4$l2g_gene)) # 27%
suzuki_st4 %>%
  select(index_snv, pos_id_b38, l2g_gene) -> st4_out
write_tsv(st4_out, 'data/t2d/suzuki_2023_st4_mapped.tsv', na='')

pp_dia_raw = read.xlsx('data_received/type_ii_03292023.xlsx')

omim_t2d = read_tsv('data/t2d/omim_t2d.tsv', col_types=cols())

vujkovic_t5 = read_tsv('data/t2d/vujkovic_supp_t5.tsv', col_types=cols())
vujkovic_t13 = suppressMessages(read_tsv('data/t2d/vujkovic_supp_t13.tsv', col_types=cols()))
vujkovic_t15 = read_tsv('data/t2d/vujkovic_supp_t15.tsv', col_types=cols())

first_novel = min(which(vujkovic_t5$description %in% c('1. novel transethnic SNP in at least one ancestral group','2. novel T2D SNP in transethnic meta only')))  
first_repl = which(vujkovic_t5$description=='3. established T2D variant')
vujkovic_t5$novel = 1:nrow(vujkovic_t5) %in% first_novel:(first_repl-1)
vujkovic_t5$gene = vujkovic_t5$nearestgene

vujkovic_t13$novel = vujkovic_t13$codingnovel==1

vujkovic_t15$best_coloc_snp_pp_h4 = suppressWarnings(as.numeric(vujkovic_t15$best_coloc_snp_pp_h4))
vujkovic_t15$pp_h4_abf = suppressWarnings(as.numeric(vujkovic_t15$pp_h4_abf))
vujkovic_t15$p = suppressWarnings(as.numeric(vujkovic_t15$p))

vujkovic_t15 %>%
  filter(pp_h4_abf >= 0.9) %>%
  filter(best_coloc_snp_pp_h4 >= 0.9) %>%
  filter(p < 5e-8) %>%
  distinct(gene) %>%
  arrange(gene) -> eqtl

eqtl$novelty = ''
eqtl$novelty[!(eqtl$gene %in% vujkovic_t5$gene[!vujkovic_t5$novel]) & !(eqtl$gene %in% vujkovic_t13$gene[!vujkovic_t13$novel])] = 'novel'
eqtl$novelty[((eqtl$gene %in% vujkovic_t5$gene[!vujkovic_t5$novel]) | (eqtl$gene %in% vujkovic_t13$gene[!vujkovic_t13$novel]))] = 'estab'

eqtl$novel = eqtl$novelty=='novel'

vujkovic_t5 %>%
  filter(gene != '-') %>%
  select(gene, novel) %>%
  mutate(source='Vujkovic 2020') -> g1
vujkovic_t13 %>%
  filter(gene != '-') %>%
  select(gene, novel) %>%
  mutate(source='Vujkovic 2020') -> g2
eqtl %>%
  select(gene, novel) %>%
  mutate(source='Vujkovic 2020')  -> g3
rbind(g1, g2, g3) -> vujkovic_genes
st4_out %>%
  select(gene=l2g_gene) %>%
  filter(!is.na(gene)) %>%
  mutate(novel = !(gene %in% vujkovic_genes$gene[!vujkovic_genes$novel])) %>%
  mutate(source='Suzuki 2023') -> suzuki_genes

t2d_gwas = rbind(vujkovic_genes, suzuki_genes) 
# now if multiple entries, and any are "established", set all to established
t2d_gwas$novel[t2d_gwas$gene %in% t2d_gwas$gene[!t2d_gwas$novel]] = FALSE
t2d_gwas %>%
  arrange(novel, desc(source)) %>%
  mutate(dup = duplicated(gene)) %>%
  filter(!dup) %>%
  mutate(source = case_when(novel ~ source,
                            !novel ~ 'Established')) %>%
  select(-dup) -> t2d_gwas

write_tsv(t2d_gwas, 'data/t2d/t2d_gwas_vujkovic_suzuki.tsv', na='')

curated = read_tsv('data/t2d/drug_gene_curation_results.tsv', col_types=cols())
pp_dia_raw$remap = curated$remap[match(pp_dia_raw$gene, curated$gene)]
pp_dia_raw$gene[!is.na(pp_dia_raw$remap)] = pp_dia_raw$remap[!is.na(pp_dia_raw$remap)]

pp_dia_raw %>%
  filter(!is.na(gene)) %>%
  left_join(t2d_gwas, by='gene') %>%
  mutate(gensup_rank = case_when(gene %in% omim_t2d$gene ~ 1, 
                                 !novel ~ 2,
                                 novel ~ 3,
                                 TRUE ~ 4)) %>% # where a drug has multiple targets, take the one with the "earliest" genetic support
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
