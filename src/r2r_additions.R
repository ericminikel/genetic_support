
mesh_best_names = read_tsv('data/mesh_best_names.tsv.gz', col_types=cols())

# all synonyms for all MeSH terms
all_vocab = read_tsv('data/mesh_all_vocab.tsv.gz', col_types=cols())
# prepare to match them to Nelson 2015
all_vocab %>%
  mutate(labeltext = tolower(labeltext)) %>%
  group_by(labeltext) %>%
  slice(1) %>%
  ungroup() -> vocab_match
# preferred terms for all MeSH IDs
best_names = read_tsv('data/mesh_best_names.tsv.gz', col_types=cols())
sim = read_tsv('data/sim.tsv.gz', col_types=cols())

# Nelson 2015 Pharmaprojects snapshot
n15p = read_tsv('data/nelson_2015_supplementary_dataset_3.tsv', col_types=cols()) %>%
  clean_names() %>%
  filter(!(phase_latest %in% c('Discontinued','No Development Reported'))) %>%
  mutate(ccat = gsub(' Clinical Trial','',phase_latest)) %>%
  mutate(mesh = tolower(msh))
n15p$mesh_id = vocab_match$id[match(n15p$msh, vocab_match$labeltext)]
n15p$mesh_term = mesh_best_names$labeltext[match(n15p$mesh_id, mesh_best_names$id)]
n15p %>%
  select(gene, mesh_id, mesh_term, ccat) %>%
  inner_join(meta_ccat, by=c('ccat'='cat')) %>%
  rename(ccatnum = num) -> n15p
# mean(!is.na(n15p$mesh_id)) # 100%
# Nelson 2015 genetics snapshot
n15g = read_tsv('data/nelson_2015_supplementary_dataset_2.tsv', col_types=cols()) %>%
  clean_names() %>%
  mutate(mesh = tolower(msh))
n15g$mesh_id = vocab_match$id[match(n15g$msh, vocab_match$labeltext)]
n15g$mesh_term = mesh_best_names$labeltext[match(n15g$mesh_id, mesh_best_names$id)]
n15g %>%
  filter(!is.na(mesh_id)) %>%
  select(gene, mesh_id, mesh_term, source) -> n15g
# mean(!is.na(n15g$mesh_id)) # 99.4%

sim_temp = sim %>%
  filter(meshcode_a %in% n15p$mesh_id) %>%
  filter(meshcode_b %in% n15g$mesh_id)

n15p %>%
  left_join(n15g, by='gene', suffix=c('_indication','_association')) %>%
  left_join(sim_temp, by=c('mesh_id_indication'='meshcode_a', 'mesh_id_association'='meshcode_b')) %>%
  mutate(comb_norm = replace_na(comb_norm, 0)) %>%
  group_by(gene, mesh_id_indication, mesh_term_indication, ccat, ccatnum) %>%
  arrange(desc(comb_norm)) %>%
  slice(1) %>%
  mutate(gensup = comb_norm >= 0.8) -> p13_g13

p13_g13 %>%
  group_by(ccatnum, ccat) %>%
  summarize(.groups='keep',
            supported = sum(gensup),
            total = n(),
            binom_obj = binom.confint(x=sum(gensup), n=n(), method='wilson')) %>%
  mutate(pg_mean = binom_obj$mean, pg_l95 = binom_obj$lower, pg_u95 = binom_obj$upper) %>%
  mutate(drug_data = 2013, genetic_data = '2013') %>%
  select(drug_data, genetic_data, ccatnum, ccat, pg_mean, pg_l95, pg_u95, supported, total) -> p13_g13_smry

assoc_temp = assoc %>%
  filter(source == 'OTG' & l2g_share >= 0.5 | source == 'PICCOLO' | source == 'OMIM' | source=='Genebass') %>%
  select(gene, mesh_id, mesh_term, source)
sim_temp = sim %>%
  filter(meshcode_a %in% n15p$mesh_id) %>%
  filter(meshcode_b %in% assoc_temp$mesh_id)
n15p %>%
  left_join(assoc_temp, by='gene', suffix=c('_indication','_association')) %>%
  left_join(sim_temp, by=c('mesh_id_indication'='meshcode_a', 'mesh_id_association'='meshcode_b')) %>%
  mutate(comb_norm = replace_na(comb_norm, 0)) %>%
  group_by(gene, mesh_id_indication, mesh_term_indication, ccat, ccatnum) %>%
  arrange(desc(comb_norm)) %>%
  slice(1) %>%
  mutate(gensup = comb_norm >= 0.8) -> p13_g23

p13_g23 %>%
  group_by(ccatnum, ccat) %>%
  summarize(.groups='keep',
            supported = sum(gensup),
            total = n(),
            binom_obj = binom.confint(x=sum(gensup), n=n(), method='wilson')) %>%
  mutate(pg_mean = binom_obj$mean, pg_l95 = binom_obj$lower, pg_u95 = binom_obj$upper) %>%
  mutate(drug_data = 2013, genetic_data = '2023') %>%
  select(drug_data, genetic_data, ccatnum, ccat, pg_mean, pg_l95, pg_u95, supported, total) -> p13_g23_smry


pp_temp = pp %>%
  rename(mesh_id = indication_mesh_id, mesh_term = indication_mesh_term) %>%
  select(gene, mesh_id, mesh_term, ccat, ccatnum)
sim_temp = sim %>%
  filter(meshcode_a %in% pp_temp$mesh_id) %>%
  filter(meshcode_b %in% n15g$mesh_id)

pp_temp %>%
  left_join(n15g, by='gene', suffix=c('_indication','_association')) %>%
  left_join(sim_temp, by=c('mesh_id_indication'='meshcode_a', 'mesh_id_association'='meshcode_b')) %>%
  mutate(comb_norm = replace_na(comb_norm, 0)) %>%
  group_by(gene, mesh_id_indication, mesh_term_indication, ccat, ccatnum) %>%
  arrange(desc(comb_norm)) %>%
  slice(1) %>%
  mutate(gensup = comb_norm >= 0.8) -> p23_g13

p23_g13 %>%
  group_by(ccatnum, ccat) %>%
  summarize(.groups='keep',
            supported = sum(gensup),
            total = n(),
            binom_obj = binom.confint(x=sum(gensup), n=n(), method='wilson')) %>%
  mutate(pg_mean = binom_obj$mean, pg_l95 = binom_obj$lower, pg_u95 = binom_obj$upper) %>%
  mutate(drug_data = 2023, genetic_data = '2013') %>%
  select(drug_data, genetic_data, ccatnum, ccat, pg_mean, pg_l95, pg_u95, supported, total) -> p23_g13_smry


sim_temp = sim %>%
  filter(meshcode_a %in% pp_temp$mesh_id) %>%
  filter(meshcode_b %in% assoc_temp$mesh_id)

pp_temp %>%
  left_join(assoc_temp, by='gene', suffix=c('_indication','_association')) %>%
  left_join(sim_temp, by=c('mesh_id_indication'='meshcode_a', 'mesh_id_association'='meshcode_b')) %>%
  mutate(comb_norm = replace_na(comb_norm, 0)) %>%
  group_by(gene, mesh_id_indication, mesh_term_indication, ccat, ccatnum) %>%
  arrange(desc(comb_norm)) %>%
  slice(1) %>%
  mutate(gensup = comb_norm >= 0.8) -> p23_g23

p23_g23 %>%
  group_by(ccatnum, ccat) %>%
  summarize(.groups='keep',
            supported = sum(gensup),
            total = n(),
            binom_obj = binom.confint(x=sum(gensup), n=n(), method='wilson')) %>%
  mutate(pg_mean = binom_obj$mean, pg_l95 = binom_obj$lower, pg_u95 = binom_obj$upper) %>%
  mutate(drug_data = 2023, genetic_data = '2023') %>%
  select(drug_data, genetic_data, ccatnum, ccat, pg_mean, pg_l95, pg_u95, supported, total) -> p23_g23_smry


# the above comparisons do not filter for genetic insight (since that is a moving target in 2013 vs. 2023)
# note that for pipeline_best though the default is to require insight, and advancement_forest always requires insight
# so to compare apples to apples, let's do the OTG pre/post 2013 comparison manually rather than using those functions

assoc_temp = assoc %>%
  filter(source == 'OTG' & l2g_share >= 0.5 & year <= 2013) %>%
  select(gene, mesh_id, mesh_term, source)
sim_temp = sim %>%
  filter(meshcode_a %in% pp_temp$mesh_id) %>%
  filter(meshcode_b %in% assoc_temp$mesh_id)

pp_temp %>%
  left_join(assoc_temp, by='gene', suffix=c('_indication','_association')) %>%
  left_join(sim_temp, by=c('mesh_id_indication'='meshcode_a', 'mesh_id_association'='meshcode_b')) %>%
  mutate(comb_norm = replace_na(comb_norm, 0)) %>%
  group_by(gene, mesh_id_indication, mesh_term_indication, ccat, ccatnum) %>%
  arrange(desc(comb_norm)) %>%
  slice(1) %>%
  mutate(gensup = comb_norm >= 0.8) -> p23_otgpre2013

p23_otgpre2013 %>%
  group_by(ccatnum, ccat) %>%
  summarize(.groups='keep',
            supported = sum(gensup),
            total = n(),
            binom_obj = binom.confint(x=sum(gensup), n=n(), method='wilson')) %>%
  mutate(pg_mean = binom_obj$mean, pg_l95 = binom_obj$lower, pg_u95 = binom_obj$upper) %>%
  mutate(drug_data = 2023, genetic_data = 'OTG through 2013') %>%
  select(drug_data, genetic_data, ccatnum, ccat, pg_mean, pg_l95, pg_u95, supported, total) -> p23_otgpre2013_smry


assoc_temp = assoc %>%
  filter(source == 'OTG' & l2g_share >= 0.5) %>%
  select(gene, mesh_id, mesh_term, source)
sim_temp = sim %>%
  filter(meshcode_a %in% pp_temp$mesh_id) %>%
  filter(meshcode_b %in% assoc_temp$mesh_id)

pp_temp %>%
  left_join(assoc_temp, by='gene', suffix=c('_indication','_association')) %>%
  left_join(sim_temp, by=c('mesh_id_indication'='meshcode_a', 'mesh_id_association'='meshcode_b')) %>%
  mutate(comb_norm = replace_na(comb_norm, 0)) %>%
  group_by(gene, mesh_id_indication, mesh_term_indication, ccat, ccatnum) %>%
  arrange(desc(comb_norm)) %>%
  slice(1) %>%
  mutate(gensup = comb_norm >= 0.8) -> p23_otgalltime

p23_otgalltime %>%
  group_by(ccatnum, ccat) %>%
  summarize(.groups='keep',
            supported = sum(gensup),
            total = n(),
            binom_obj = binom.confint(x=sum(gensup), n=n(), method='wilson')) %>%
  mutate(pg_mean = binom_obj$mean, pg_l95 = binom_obj$lower, pg_u95 = binom_obj$upper) %>%
  mutate(drug_data = 2023, genetic_data = 'OTG all time') %>%
  select(drug_data, genetic_data, ccatnum, ccat, pg_mean, pg_l95, pg_u95, supported, total) -> p23_otgalltime_smry


rbind(p13_g13_smry, p13_g23_smry, p23_g13_smry, p23_g23_smry, p23_otgpre2013_smry, p23_otgalltime_smry) -> pg_time_comparison
# write_supp_table(pg_time_comparison, 'P(G) by phase using 2013 vs. 2023 drug and genetics datasets.')


pg_time_comparison %>%
  mutate(y = 6-ccatnum) %>%
  rename(mean=pg_mean, l95=pg_l95, u95=pg_u95, numerator=supported, denominator=total, label=ccat) -> pg_forest

resx=300
png('display_items/figure_s_new1.png',width=6.5*resx,height=4.5*resx,res=resx)
par(mfrow=c(3,2))
panel = 1
plot_forest(pg_forest %>% filter(drug_data==2013 & genetic_data=='2013'), title='2013 drug pipeline\n2013 genetics', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5)
panel = panel + 1
plot_forest(pg_forest %>% filter(drug_data==2013 & genetic_data=='2023'), title='2013 drug pipeline\n2023 genetics', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5)
panel = panel + 1
plot_forest(pg_forest %>% filter(drug_data==2023 & genetic_data=='2013'), title='2023 drug pipeline\n2013 genetics', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5)
panel = panel + 1
plot_forest(pg_forest %>% filter(drug_data==2023 & genetic_data=='2023'), title='2023 drug pipeline\n2023 genetics', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5)
panel = panel + 1
plot_forest(pg_forest %>% filter(drug_data==2023 & genetic_data=='OTG through 2013'), title='2023 drug pipeline\nOTG only, 2005-2013', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5)
panel = panel + 1
plot_forest(pg_forest %>% filter(drug_data==2023 & genetic_data=='OTG all time'), title='2023 drug pipeline\nOTG only, all time', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5)
panel = panel + 1
dev.off()



# cross-tab of 2013 area assignments vs. 2023 area assignments
n15p = read_tsv('data/nelson_2015_supplementary_dataset_3.tsv', col_types=cols()) %>%
  clean_names() %>%
  filter(!(phase_latest %in% c('Discontinued','No Development Reported'))) %>%
  mutate(ccat = gsub(' Clinical Trial','',phase_latest)) %>%
  mutate(mesh = tolower(msh))
n15p$mesh_id = vocab_match$id[match(n15p$msh, vocab_match$labeltext)]
n15p$mesh_term = mesh_best_names$labeltext[match(n15p$mesh_id, mesh_best_names$id)]

n15p %>%
  distinct(msh_top) -> n15_areas

#hmmm. there are 25 msh_top values which do not exactly match the ones used int he Nelson 2015 figures.
# there must exist another layer of mapping.

# pull out the therapeutic areas (called "categories" then) from Figure 3A in Nelson 2015
n15a = read_tsv('data/nelson_2015_table_s6.tsv', col_types=cols()) %>%
  clean_names() %>%
  distinct(msh_ind, category) 

cats = read_tsv('data/nelson_2015_categories.tsv', col_types=cols()) %>%
  mutate(y = max(ord) - ord + 1)

read_tsv('data/areas.tsv', col_types=cols()) %>%
  select(topl, area, color) %>%
  mutate(x=row_number()) -> current_areas

n15a %>%
  inner_join(vocab_match, by=c('msh_ind'='labeltext')) %>%
  inner_join(indic_topl_match, by=c('id'='indication_mesh_id')) %>%
  inner_join(current_areas, by='topl') %>%
  inner_join(cats, by='category') %>%
  group_by(y, category, x, area, color) %>%
  summarize(.groups='keep', n = length(unique(id))) %>%
  ungroup() -> area_xtab

area_xtab$plot_color = alpha(area_xtab$color, area_xtab$n / max(area_xtab$n))

resx=300
png('display_items/figure_s_new2.png',width=3.25*resx,height=3.25*resx,res=resx)
par(mar=c(7,7,1,1))
xlims = range(current_areas$x) + c(-0.6, 0.6)
ylims = range(cats$y) + c(-0.6, 0.6)
boxrad = 0.5
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
abline(v=unique(c(boxrad, current_areas$x + boxrad)), lwd=0.125)
abline(h=unique(c(boxrad, cats$y + boxrad)), lwd=0.125)
rect(xleft=area_xtab$x - boxrad, xright=area_xtab$x + boxrad, ybottom = area_xtab$y - boxrad, ytop = area_xtab$y + boxrad,
     col = area_xtab$plot_color, border='#000000', lwd=0.5)
text(x=area_xtab$x, y=area_xtab$y, labels=area_xtab$n, cex=0.5)
mtext(side=1, at=current_areas$x, text=current_areas$area, col=current_areas$color, las=2, cex=0.6)
mtext(side=2, at=cats$y, text=cats$category, col='#4D4D4D', las=2, cex=0.6)
dev.off()



# "Add a new supp figure of P(G) vs. Phase (like Fig 1A) for each of the 20 therapy areas. One panel per therapy area."
resx=300
png('display_items/figure_s_new3.png',width=6.5*resx,height=8.0*resx,res=resx)
panel = 1
par(mfrow = c(6,3))
for (i in 1:nrow(areas_all)) {
  
  combined_ti_area = subset_by_area(combined_ti, topl=areas_all$topl[i], filter=areas_all$filter[i])
  area_forest = advancement_forest(combined_ti_area)
  plot_forest(area_forest, xlims=c(0,.60), xlab='P(G) vs. phase', col=areas_all$color[i], title='', mar=c(3,6,2.5,5))
  mtext(side=3, adj=0, at=-0.3, line=0.25, col='#000000', text=areas_all$area[i], cex=.85)
  mtext(letters[panel], side=3, cex=1.4, at=-0.5, line = 0.25)
  panel = panel + 1
  
  area_forest$area = areas_all$area[i]
  if (i == 1) {
    pg_out = area_forest
  } else {
    pg_out = rbind(pg_out, area_forest)
  }
}
dev.off()

pg_out %>%
  select(area, phase=label, phasenum=num, supported=numerator, total=denominator, pg_mean=mean, pg_l95 = l95, pg_u95 = u95) -> pg_out

write_supp_table(pg_out, "Proportion of target-indication pairs with genetic support by phase and therapy area, combined mode (both historical and active programs).")


##### stuff for rev 3
assoc %>%
  filter(source=='PICCOLO') %>%
  mutate(qtl_source = case_when(grepl('gtex',tolower(extra_info)) ~ 'GTEx',
                                TRUE ~ 'other')) %>%
  group_by(qtl_source) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup()

assoc %>%
  filter(source=='OTG') %>%
  mutate(subsource = gsub('[0-9_].*','',gsub('https://genetics.opentargets.org/study/','',original_link))) %>%
  group_by(subsource) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup()

##### confounding between TA and categories in panel 1D
# yearfirst
for (i in 1:nrow(areas)) {
  this_topl = areas$topl[i]
  this_area = areas$area[i]
  yf_obj = subset_by_area(yearfirst_logit_data, topl=this_topl, filter='only')
  yf_minimal = yf_obj %>%
    mutate(area = this_area) %>%
    mutate(topl = this_topl) %>%
    select(topl, area, uid, gene, indication_mesh_id, indication_mesh_term, assoc_year)
  if (i == 1) {
    area_x_year = yf_minimal
  } else {
    area_x_year = rbind(yf_minimal, area_x_year)
  }
}
# Kruskal Wallis test
area_year_kw = kruskal.test(assoc_year ~ area, data=area_x_year)
area_year_kw_p = area_year_kw$p.value
# also try a discrete model more analogous to what 1D shows visually
area_x_year$yearbin = floor((area_x_year$assoc_year - 2007)/4)*4+2007
area_year_ctable = table(area_x_year[,c('yearbin','area')])
area_year_chisq_p = suppressWarnings(chisq.test(area_year_ctable))$p.value

write(paste("Confounding between therapy area and year of discovery among OTG GWAS Catalog-supported T-I pairs: KW test P = ",
            formatC(area_year_kw_p, format='e', digits=1),
            ', Chi Square test P = ',formatC(area_year_chisq_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)

# genecount
for (i in 1:nrow(areas)) {
  this_topl = areas$topl[i]
  this_area = areas$area[i]
  gc_obj = subset_by_area(gc_logit_data, topl=this_topl, filter='only')
  gc_minimal = gc_obj %>%
    mutate(area = this_area) %>%
    mutate(topl = this_topl) %>%
    select(topl, area, uid, gene, indication_mesh_id, indication_mesh_term, gene_count)
  if (i == 1) {
    area_x_gc = gc_minimal
  } else {
    area_x_gc = rbind(gc_minimal, area_x_gc)
  }
}
# KW
area_gc_kw = kruskal.test(gene_count ~ area, data=area_x_gc)
area_gc_kw_p = area_gc_kw$p.value
# also try a discrete model more analogous to what 1D shows visually
#10^floor(log10(area_x_gc$gene_count))
area_x_gc$gcbin = case_when(area_x_gc$gene_count == 1 ~ 1,
                            area_x_gc$gene_count >= 2 & area_x_gc$gene_count <= 9 ~ 2,
                            area_x_gc$gene_count >= 10 & area_x_gc$gene_count <= 99 ~ 10,
                            area_x_gc$gene_count >= 100 & area_x_gc$gene_count <= 999 ~ 100,
                            area_x_gc$gene_count >= 1000 ~ 1000)
area_gc_ctable = table(area_x_gc[,c('gcbin','area')])
area_gc_chisq_p = suppressWarnings(chisq.test(area_gc_ctable)$p.value)

write(paste("Confounding between therapy area and gene count among OTG GWAS Catalog-supported T-I pairs: KW test P = ",
            formatC(area_gc_kw_p, format='e', digits=1),
            ', Chi Square test P = ',formatC(area_gc_chisq_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)


# beta
for (i in 1:nrow(areas)) {
  this_topl = areas$topl[i]
  this_area = areas$area[i]
  beta_obj = subset_by_area(beta_logit_data, topl=this_topl, filter='only')
  beta_minimal = beta_obj %>%
    mutate(area = this_area) %>%
    mutate(topl = this_topl) %>%
    select(topl, area, uid, gene, indication_mesh_id, indication_mesh_term, abs_beta)
  if (i == 1) {
    area_x_beta = beta_minimal
  } else {
    area_x_beta = rbind(beta_minimal, area_x_beta)
  }
}
# KW
area_beta_kw = kruskal.test(abs_beta ~ area, data=area_x_beta)
area_beta_kw_p = area_beta_kw$p.value
# also discrete model like 1D
crossing(area_x_beta, (beta_rrs %>% filter(label!='All with beta'))) %>%
           filter(abs_beta >= min_beta & abs_beta <= max_beta) -> area_x_beta_binned
area_beta_ctable = table(area_x_beta_binned[,c('label','area')])
area_beta_chisq_p = suppressWarnings(chisq.test(area_beta_ctable)$p.value)

write(paste("Confounding between therapy area and beta among OTG GWAS Catalog-supported T-I pairs: KW test P = ",
            formatC(area_beta_kw_p, format='e', digits=1),
            ', Chi Square test P = ',formatC(area_beta_chisq_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)


# or
for (i in 1:nrow(areas)) {
  this_topl = areas$topl[i]
  this_area = areas$area[i]
  or_obj = subset_by_area(or_logit_data, topl=this_topl, filter='only')
  or_minimal = or_obj %>%
    mutate(area = this_area) %>%
    mutate(topl = this_topl) %>%
    select(topl, area, uid, gene, indication_mesh_id, indication_mesh_term, abs_or)
  if (i == 1) {
    area_x_or = or_minimal
  } else {
    area_x_or = rbind(or_minimal, area_x_or)
  }
}
# KW
area_or_kw = kruskal.test(abs_or ~ area, data=area_x_or)
area_or_kw_p = area_or_kw$p.value
# also discrete model like 1D
crossing(area_x_or, (or_rrs %>% filter(label!='All with OR'))) %>%
  filter(abs_or >= min_or & abs_or <= max_or) -> area_x_or_binned
area_or_ctable = table(area_x_or_binned[,c('label','area')])
area_or_chisq_p = suppressWarnings(chisq.test(area_or_ctable)$p.value)

write(paste("Confounding between therapy area and OR among OTG GWAS Catalog-supported T-I pairs: KW test P = ",
            formatC(area_or_kw_p, format='e', digits=1),
            ', Chi Square test P = ',formatC(area_or_chisq_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)



# MAF
for (i in 1:nrow(areas)) {
  this_topl = areas$topl[i]
  this_area = areas$area[i]
  maf_obj = subset_by_area(maf_logit_data, topl=this_topl, filter='only')
  maf_minimal = maf_obj %>%
    mutate(area = this_area) %>%
    mutate(topl = this_topl) %>%
    select(topl, area, uid, gene, indication_mesh_id, indication_mesh_term, lead_maf)
  if (i == 1) {
    area_x_maf = maf_minimal
  } else {
    area_x_maf = rbind(maf_minimal, area_x_maf)
  }
}
# KW maf continuous model
area_maf_kw = kruskal.test(lead_maf ~ area, data=area_x_maf)
area_maf_kw_p = area_maf_kw$p.value
# also discrete model like 1D
crossing(area_x_maf, (maf_rrs %>% filter(label!='All with MAF'))) %>%
  filter(lead_maf >= min_maf & lead_maf <= max_maf) -> area_x_maf_binned
area_maf_ctable = table(area_x_maf_binned[,c('label','area')])
area_maf_chisq_p = suppressWarnings(chisq.test(area_maf_ctable)$p.value)


write(paste("Confounding between therapy area and MAF among OTG GWAS Catalog-supported T-I pairs: KW test P = ",
            formatC(area_maf_kw_p, format='e', digits=1),
            ', Chi Square test P = ',formatC(area_maf_chisq_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)




### assess confounding between TA and OTG source
for (i in 1:nrow(areas)) {
  this_topl = areas$topl[i]
  this_area = areas$area[i]
  gwascat_obj = subset_by_area(hist_ti_gwascat, topl=this_topl, filter='only')
  ukbb_obj = subset_by_area(hist_ti_ukbb, topl=this_topl, filter='only')
  finngen_obj = subset_by_area(hist_ti_finngen, topl=this_topl, filter='only')
  rbind(gwascat_obj, ukbb_obj, finngen_obj) %>%
    filter(target_status == 'genetically supported target' & !is.na(indication_mesh_id) & cat %in% c('Launched',active_clinical$cat)) %>%
    mutate(area = this_area) %>%
    mutate(topl = this_topl) %>%
    select(topl, area, uid, gene, indication_mesh_id, indication_mesh_term, gwas_source) -> subcat_minimal
  
  if (i == 1) {
    area_x_subcat = subcat_minimal
  } else {
    area_x_subcat = rbind(subcat_minimal, area_x_subcat)
  }
}

area_subcat_ctable = table(area_x_subcat[,c('gwas_source','area')])
area_subcat_chisq_p = suppressWarnings(chisq.test(area_subcat_ctable)$p.value)



write(paste("Confounding between therapy area and GWAS source: ",
            'Chi Square test P = ',formatC(area_subcat_chisq_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)


resx=300
png(paste0(output_path,'/figure-s6.png'),width=6.5*resx,height=4*resx,res=resx)

layout_matrix = matrix(1:7, byrow=T, nrow=1)
layout(layout_matrix)

if (is.null(areas$y)) {
  areas$y = nrow(areas):1
}
line_color = '#B9B9B9'
xlims = c(0,5)
ylims = range(areas$y, na.rm=T) + c(-0.5, 0.5)
par(mar=c(3,1,3,0))
plot(NA, NA, ylim=ylims, xlim=xlims, axes=F, ann=F, xaxs='i', yaxs='i')
mtext(side=4, at=areas$y, line=0, adj=1, las=2, text=areas$area, cex=.7, col=areas$color)
panel = 1
par(mar=c(3,0.25,3,0.75))
xlims = c(2007, 2022) + c(-0.5, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, at=2007:2022, tck=-0.025, labels=NA)
axis(side=1, at=seq(2007, 2022, 4), tck=-0.05, labels=NA)
axis(side=1, at=seq(2007, 2022, 4), line=-0.75, lwd=0, cex.axis=0.7)
abline(v=min(xlims))
area_x_year %>%
  inner_join(areas, by='topl') %>%
  select(x=assoc_year, y, color) -> toplot
set.seed(1)
points(jitter(toplot$x,amount=.25), jitter(toplot$y,amount=.25), col=alpha(toplot$color, .1), pch=20)
toplot %>%
  group_by(y, color) %>%
  summarize(.groups='keep', 
            median_x = median(x),
            q25 = quantile(x,.25),
            q75 = quantile(x,.75)) %>%
  ungroup() -> tosmry
barwidth = .33
segments(x0=tosmry$median_x, y0=tosmry$y - barwidth, y1=tosmry$y + barwidth, lwd=2, lend=1, col=tosmry$color)
rect(xleft=tosmry$q25, xright=tosmry$q75, ybottom=tosmry$y - barwidth, ytop=tosmry$y + barwidth, lwd=1.5, border=tosmry$color, col=NULL)
mtext(side=1, cex=0.8, line=1.3, text='year')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


xlims = c(1,2000)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
xats = rep(1:9, 4) * 10^rep(0:3, each=9)
xbigs = 10^(0:3)
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, line=-0.75, lwd=0, cex.axis=0.7)
abline(v=min(xlims))
area_x_gc %>%
  inner_join(areas, by='topl') %>%
  select(x=gene_count, y, color) -> toplot
set.seed(1)
points(jitter(toplot$x,amount=.25), jitter(toplot$y,amount=.25), col=alpha(toplot$color, .1), pch=20)
toplot %>%
  group_by(y, color) %>%
  summarize(.groups='keep', 
            median_x = median(x),
            q25 = quantile(x,.25),
            q75 = quantile(x,.75)) %>%
  ungroup() -> tosmry
barwidth = .33
segments(x0=tosmry$median_x, y0=tosmry$y - barwidth, y1=tosmry$y + barwidth, lwd=2, lend=1, col=tosmry$color)
rect(xleft=tosmry$q25, xright=tosmry$q75, ybottom=tosmry$y - barwidth, ytop=tosmry$y + barwidth, lwd=1.5, border=tosmry$color, col=NULL)
mtext(side=1, cex=0.8, line=1.3, text='gene count')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


xlims = c(0.001, 50)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
xats = rep(1:9, 5) * 10^rep(-3:1, each=9)
xbigs = 10^(-3:1)
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, labels=formatC(xbigs, format='g'), line=-0.75, lwd=0, cex.axis=0.7)
abline(v=min(xlims))
area_x_beta %>%
  filter(abs_beta > 0) %>%
  inner_join(areas, by='topl') %>%
  select(x=abs_beta, y, color) -> toplot
set.seed(1)
points(jitter(toplot$x,amount=.25), jitter(toplot$y,amount=.25), col=alpha(toplot$color, .1), pch=20)
toplot %>%
  group_by(y, color) %>%
  summarize(.groups='keep', 
            median_x = median(x),
            q25 = quantile(x,.25),
            q75 = quantile(x,.75)) %>%
  ungroup() -> tosmry
barwidth = .33
segments(x0=tosmry$median_x, y0=tosmry$y - barwidth, y1=tosmry$y + barwidth, lwd=2, lend=1, col=tosmry$color)
rect(xleft=tosmry$q25, xright=tosmry$q75, ybottom=tosmry$y - barwidth, ytop=tosmry$y + barwidth, lwd=1.5, border=tosmry$color, col=NULL)
mtext(side=1, cex=0.8, line=1.3, text='beta')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


xlims = c(0.001, 50)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
xats = rep(1:9, 5) * 10^rep(-3:1, each=9)
xbigs  = c(10^(-3:-1), 1, 9)
xbigslabs = c(1+10^(-3:-1), 2, 10)
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, labels=formatC(xbigslabs, format='g'), line=-0.75, lwd=0, cex.axis=0.7)
abline(v=min(xlims))
area_x_or %>%
  inner_join(areas, by='topl') %>%
  mutate(x = abs_or-1) %>%
  select(x, y, color) -> toplot
set.seed(1)
points(jitter(toplot$x,amount=.25), jitter(toplot$y,amount=.25), col=alpha(toplot$color, .1), pch=20)
toplot %>%
  group_by(y, color) %>%
  summarize(.groups='keep', 
            median_x = median(x),
            q25 = quantile(x,.25),
            q75 = quantile(x,.75)) %>%
  ungroup() -> tosmry
barwidth = .33
segments(x0=tosmry$median_x, y0=tosmry$y - barwidth, y1=tosmry$y + barwidth, lwd=2, lend=1, col=tosmry$color)
rect(xleft=tosmry$q25, xright=tosmry$q75, ybottom=tosmry$y - barwidth, ytop=tosmry$y + barwidth, lwd=1.5, border=tosmry$color, col=NULL)
mtext(side=1, cex=0.8, line=1.3, text='OR')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


xlims = c(0.001, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
xats = rep(1:9, 4) * 10^rep(-3:0, each=9)
xbigs = c(10^(-3:-1),.5)
xbigslabs = c('0.1%','1%','10%','50%')
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, labels=formatC(xbigslabs, format='g'), line=-0.75, lwd=0, cex.axis=0.7)
abline(v=min(xlims))
area_x_maf %>%
  filter(lead_maf > 0) %>%
  inner_join(areas, by='topl') %>%
  mutate(x = lead_maf) %>%
  select(x, y, color) -> toplot
set.seed(1)
points(jitter(toplot$x,amount=.25), jitter(toplot$y,amount=.25), col=alpha(toplot$color, .1), pch=20)
toplot %>%
  group_by(y, color) %>%
  summarize(.groups='keep', 
            median_x = median(x),
            q25 = quantile(x,.25),
            q75 = quantile(x,.75)) %>%
  ungroup() -> tosmry
barwidth = .33
segments(x0=tosmry$median_x, y0=tosmry$y - barwidth, y1=tosmry$y + barwidth, lwd=2, lend=1, col=tosmry$color)
rect(xleft=tosmry$q25, xright=tosmry$q75, ybottom=tosmry$y - barwidth, ytop=tosmry$y + barwidth, lwd=1.5, border=tosmry$color, col=NULL)
mtext(side=1, cex=0.8, line=1.3, text='MAF')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1



area_x_subcat %>%
  group_by(topl, area, gwas_source) %>%
  summarize(.groups='keep', n_ti=n()) %>%
  ungroup() %>%
  group_by(gwas_source) %>%
  mutate(area_proportion = n_ti/sum(n_ti)) %>%
  ungroup() -> area_x_subcat_proportions

subcat_meta = tibble(gwas_source=c('GWAS Catalog','Neale UKBB','FinnGen'),
                     x=1:3)

area_x_subcat_proportions %>%
  inner_join(select(areas, topl, area, color, y), by=c('topl','area')) %>%
  inner_join(subcat_meta, by='gwas_source') %>%
  arrange(y) %>%
  group_by(x, gwas_source) %>%
  mutate(cume_prop = cumsum(area_proportion)) %>%
  ungroup() %>%
  arrange(x, desc(y)) -> area_subcat_plotdata

par(mar=c(7,3,3,0.75))
xlims = c(0.5, 3.5)
ylims = c(0, 1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, at=subcat_meta$x, text=subcat_meta$gwas_source, cex=0.65, line=0.25, las=2)
axis(side=2, at=0:4/4, labels=NA)
axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2, lwd=0, line=-0.25, cex.axis=0.7)
barwidth=0.33
rect(xleft=area_subcat_plotdata$x - barwidth,
     xright=area_subcat_plotdata$x + barwidth,
     ybottom=rep(0,length(area_subcat_plotdata)),
     ytop=area_subcat_plotdata$cume_prop,
     col=area_subcat_plotdata$color, border=NA)

dev.off()





combined_ti_gwas_sans_omim    = pipeline_best(merge2, phase='combined', basis='ti', associations=c('PICCOLO','OTG','Genebass'), lacking=c('OMIM'), verbose=F)
combined_ti_omim_sans_gwas    = pipeline_best(merge2, phase='combined', basis='ti', associations=c('OMIM'), lacking=c('PICCOLO','OTG','Genebass'), verbose=F)
combined_ti_gwas_andalso_omim = pipeline_best(merge2, phase='combined', basis='ti', associations=c('PICCOLO','OTG','Genebass'), andalso=c('OMIM'), verbose=F)
combined_ti_omim_andalso_gwas = pipeline_best(merge2, phase='combined', basis='ti', associations=c('OMIM'), andalso=c('PICCOLO','OTG','Genebass'), verbose=F)


advancement_rr(combined_ti_gwas_sans_omim)
advancement_rr(combined_ti_omim_sans_gwas)
advancement_rr(combined_ti_omim)
advancement_rr(combined_ti_gwas_andalso_omim) # 3.93
advancement_rr(combined_ti_omim_andalso_gwas) # 3.93
# above RS numbers now agree; check that there are no T-I on which they disagree:
omim_plus_gwas_ti = unique(combined_ti_omim_andalso_gwas$ti_uid[combined_ti_omim_andalso_gwas$target_status=='genetically supported target'])
gwas_plus_omim_ti = unique(combined_ti_gwas_andalso_omim$ti_uid[combined_ti_gwas_andalso_omim$target_status=='genetically supported target'])
setdiff(omim_plus_gwas_ti, gwas_plus_omim_ti)
setdiff(gwas_plus_omim_ti, omim_plus_gwas_ti)

# experiment with historical being T-I no longer in active
pp$hcat[!is.na(pp$acat)] = NA
pp$hcatnum[!is.na(pp$acat)] = NA

merge2$hcat[!is.na(merge2$acat)] = NA
merge2$hcatnum[!is.na(merge2$acat)] = NA



combined_ti_genebass_misskat  = pipeline_best(merge2, phase='combined', basis='ti', associations='Genebass', genebass_subcat='missenseLC skat', verbose=F)
combined_ti_genebass_misburd  = pipeline_best(merge2, phase='combined', basis='ti', associations='Genebass', genebass_subcat='missenseLC burden', verbose=F)
combined_ti_genebass_lofskat  = pipeline_best(merge2, phase='combined', basis='ti', associations='Genebass', genebass_subcat='pLoF skat', verbose=F)
combined_ti_genebass_lofburd  = pipeline_best(merge2, phase='combined', basis='ti', associations='Genebass', genebass_subcat='pLoF burden', verbose=F)

advancement_rr(combined_ti_genebass_misskat)
advancement_rr(combined_ti_genebass_misburd)
advancement_rr(combined_ti_genebass_lofskat)
advancement_rr(combined_ti_genebass_lofburd)
