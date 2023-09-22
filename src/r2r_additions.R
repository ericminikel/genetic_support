
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

# examples of what is in signs & symptoms
pp %>% inner_join(indic_topl_match, by=c('indication_mesh_id','indication_mesh_term')) %>% filter(topl=='C23') %>%
  group_by(indication_mesh_id, indication_mesh_term) %>%
  summarize(.groups='keep', n=n()) %>%
  arrange(desc(n))


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




