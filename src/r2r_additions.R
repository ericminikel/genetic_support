
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
n15p$mesh_term = best_names$labeltext[match(n15p$mesh_id, best_match$id)]
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
n15g$mesh_term = best_names$labeltext[match(n15g$mesh_id, best_match$id)]
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
  mutate(drug_data = 2013, genetic_data = 2013) %>%
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
  mutate(drug_data = 2013, genetic_data = 2023) %>%
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
  mutate(drug_data = 2023, genetic_data = 2013) %>%
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
  mutate(drug_data = 2023, genetic_data = 2023) %>%
  select(drug_data, genetic_data, ccatnum, ccat, pg_mean, pg_l95, pg_u95, supported, total) -> p23_g23_smry

rbind(p13_g13_smry, p13_g23_smry, p23_g13_smry, p23_g23_smry) -> pg_time_comparison

# write_supp_table(pg_time_comparison, 'P(G) by phase using 2013 vs. 2023 drug and genetics datasets.')

pg_time_comparison %>%
  mutate(y = 6-ccatnum) %>%
  rename(mean=pg_mean, l95=pg_l95, u95=pg_u95, numerator=supported, denominator=total, label=ccat) -> pg_forest

resx=300
png('display_items/figure_s_new1.png',width=6.5*resx,height=4.5*resx,res=resx)
par(mfrow=c(2,2))
plot_forest(pg_forest %>% filter(drug_data==2013 & genetic_data==2013), title='2013 drug pipeline\n2013 genetics', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
plot_forest(pg_forest %>% filter(drug_data==2013 & genetic_data==2023), title='2013 drug pipeline\n2023 genetics', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
plot_forest(pg_forest %>% filter(drug_data==2023 & genetic_data==2013), title='2023 drug pipeline\n2013 genetics', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
plot_forest(pg_forest %>% filter(drug_data==2023 & genetic_data==2023), title='2023 drug pipeline\n2023 genetics', xlims=c(0,.12), xlab='P(G) vs. phase', col='#000000')
dev.off()

hist_ti_otg_pre2013 = pipeline_best(merge2, phase='historical', basis='ti', associations=c('OTG'), verbose=F, max_year=2013)
active_ti_otg_pre2013 = pipeline_best(merge2, phase='active', basis='ti', associations=c('OTG'), verbose=F, max_year=2013)
combined_ti_otg_pre2013 = pipeline_best(merge2, phase='combined', basis='ti', associations=c('OTG'), verbose=F, max_year=2013)

hist_ti_forest = advancement_forest(hist_ti_otg_pre2013,phase='historical')
active_ti_forest = advancement_forest(active_ti_otg_pre2013,phase='active')
combined_ti_forest = advancement_forest(combined_ti_otg_pre2013,phase='combined')

resx=300
png('display_items/r2r_pg_otg_through_2013.png',width=5.5*resx,height=3.5*resx,res=resx)
panel = 1
#### 1A - T-I pair forest
hist_col = '#F46D43'
active_col = '#74ADD1'
combined_col = '#9970AB'
hist_offset = 0.25
active_offset = 0.0
combined_offset = -0.25
plot_forest(combined_ti_forest, xlims=c(0,.15), xlab='P(G) vs. phase', col='#00000000')
mtext(side=2, at=combined_ti_forest$y, text=combined_ti_forest$label, line=0.5, las=2, cex=0.75)
points(hist_ti_forest$mean[1:4], hist_ti_forest$y[1:4] + hist_offset, pch=19, col=hist_col)
segments(x0=hist_ti_forest$l95[1:4], x1=hist_ti_forest$u95[1:4], y0=hist_ti_forest$y[1:4] + hist_offset, lwd=2, col=hist_col) # 95%CIs
points(active_ti_forest$mean[2:5], active_ti_forest$y[2:5] + active_offset, pch=19, col=active_col)
segments(x0=active_ti_forest$l95[2:5], x1=active_ti_forest$u95[2:5], y0=active_ti_forest$y[2:5] + active_offset, lwd=2, col=active_col) # 95%CIs
points(combined_ti_forest$mean[1:4], combined_ti_forest$y[1:4] + combined_offset, pch=19, col=combined_col)
segments(x0=combined_ti_forest$l95[1:4], x1=combined_ti_forest$u95[1:4], y0=combined_ti_forest$y[1:4] + combined_offset, lwd=2, col=combined_col) # 95%CIs
points(combined_ti_forest$mean[5], combined_ti_forest$y[5], pch=19, col=combined_col)
segments(x0=combined_ti_forest$l95[5], x1=combined_ti_forest$u95[5], y0=combined_ti_forest$y[5], lwd=2, col=combined_col) # 95%CIs
abline(h=0:5+.5, lwd=.5)
par(xpd=T)
legend(x=0.07, y=5.5, legend=c('historical','active','combined'), pch=15, col=c(hist_col,active_col,combined_col), text.col=c(hist_col,active_col,combined_col), cex=0.75, bg='#FFFFFF')
par(xpd=F)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1
mtext(side=3, text='OTG through 2013')
dev.off()

hist_ti_otg_all = pipeline_best(merge2, phase='historical', basis='ti', associations=c('OTG'), verbose=F)
active_ti_otg_all = pipeline_best(merge2, phase='active', basis='ti', associations=c('OTG'), verbose=F)
combined_ti_otg_all = pipeline_best(merge2, phase='combined', basis='ti', associations=c('OTG'), verbose=F)

hist_ti_forest = advancement_forest(hist_ti_otg_all,phase='historical')
active_ti_forest = advancement_forest(active_ti_otg_all,phase='active')
combined_ti_forest = advancement_forest(combined_ti_otg_all,phase='combined')


resx=300
png('display_items/r2r_pg_otg_alltime.png',width=5.5*resx,height=3.5*resx,res=resx)
panel = 1
#### 1A - T-I pair forest
hist_col = '#F46D43'
active_col = '#74ADD1'
combined_col = '#9970AB'
hist_offset = 0.25
active_offset = 0.0
combined_offset = -0.25
plot_forest(combined_ti_forest, xlims=c(0,.15), xlab='P(G) vs. phase', col='#00000000')
mtext(side=2, at=combined_ti_forest$y, text=combined_ti_forest$label, line=0.5, las=2, cex=0.75)
points(hist_ti_forest$mean[1:4], hist_ti_forest$y[1:4] + hist_offset, pch=19, col=hist_col)
segments(x0=hist_ti_forest$l95[1:4], x1=hist_ti_forest$u95[1:4], y0=hist_ti_forest$y[1:4] + hist_offset, lwd=2, col=hist_col) # 95%CIs
points(active_ti_forest$mean[2:5], active_ti_forest$y[2:5] + active_offset, pch=19, col=active_col)
segments(x0=active_ti_forest$l95[2:5], x1=active_ti_forest$u95[2:5], y0=active_ti_forest$y[2:5] + active_offset, lwd=2, col=active_col) # 95%CIs
points(combined_ti_forest$mean[1:4], combined_ti_forest$y[1:4] + combined_offset, pch=19, col=combined_col)
segments(x0=combined_ti_forest$l95[1:4], x1=combined_ti_forest$u95[1:4], y0=combined_ti_forest$y[1:4] + combined_offset, lwd=2, col=combined_col) # 95%CIs
points(combined_ti_forest$mean[5], combined_ti_forest$y[5], pch=19, col=combined_col)
segments(x0=combined_ti_forest$l95[5], x1=combined_ti_forest$u95[5], y0=combined_ti_forest$y[5], lwd=2, col=combined_col) # 95%CIs
abline(h=0:5+.5, lwd=.5)
par(xpd=T)
legend(x=0.07, y=5.5, legend=c('historical','active','combined'), pch=15, col=c(hist_col,active_col,combined_col), text.col=c(hist_col,active_col,combined_col), cex=0.75, bg='#FFFFFF')
par(xpd=F)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1
mtext(side=3, text='OTG all time')
dev.off()



pp %>% inner_join(indic_topl_match, by=c('indication_mesh_id','indication_mesh_term')) %>% filter(topl=='C23') %>%
  group_by(indication_mesh_id, indication_mesh_term) %>%
  summarize(.groups='keep', n=n()) %>%
  arrange(desc(n))
