

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

# for the Figure S6 stuff, want to be able to calculate RS based solely on
# columns ccat, ccatnum, and gensup
# example: p13_g13
adv_rr_simple = function(ptbl, alpha=0.05) {
  
  phase_map = tibble(phase=factor(c('Preclinical','I','II','III','I-Launch'),ordered=T,levels=c('Preclinical','I','II','III','I-Launch')), 
                     phorder = 0:4,
                     varname=c('succ_p_1','succ_1_2','succ_2_3','succ_3_a','succ_1_a'))
  rr_long_structure = crossing(gs=c('yes','no'),phase=c('Preclinical','I','II','III'))
  ptbl %>% 
    mutate(ti_uid = paste0(gene,'-',mesh_id_indication)) %>%
    filter(ccat %in% c('Preclinical','Phase I','Phase II','Phase III','Launched')) %>%
    mutate(succ_p_1 = case_when(ccat %in% c('Phase I','Phase II','Phase III','Launched') ~ TRUE,
                                ccat %in% c('Preclinical') ~ FALSE)) %>%
    mutate(succ_1_2 = case_when(ccat %in% c('Phase II','Phase III','Launched') ~ TRUE,
                                ccat %in% c('Phase I') ~ FALSE,
                                TRUE ~ NA)) %>%
    mutate(succ_2_3 = case_when(ccat %in% c('Phase III','Launched') ~ TRUE,
                                ccat %in% c('Phase II') ~ FALSE,
                                TRUE ~ NA)) %>%
    mutate(succ_3_a = case_when(ccat %in% c('Launched') ~ TRUE,
                                ccat %in% c('Phase III') ~ FALSE,
                                TRUE ~ NA)) %>%
    mutate(gs = ifelse(gensup,'yes','no')) %>%
    select(gs, ti_uid, ccatnum, ccat, succ_p_1, succ_1_2, succ_2_3, succ_3_a) %>%
    pivot_longer(succ_p_1:succ_3_a) %>%
    rename(success=value) %>%
    filter(!is.na(success)) %>%
    inner_join(phase_map, by=c('name'='varname')) %>%
    select(ti_uid, gs, phase, success) -> long
  
  long %>%
    group_by(gs, phase) %>%
    summarize(.groups='keep',
              x = sum(success),
              n = sum(!is.na(success))) %>%
    ungroup() %>%
    right_join(rr_long_structure, by=c('gs','phase')) %>%
    mutate(x=replace_na(x, 0),
           n=replace_na(n, 0)) %>%
    mutate(binom = binom.confint(x, n, 1-alpha, method='wilson')[,c('mean','lower','upper')]) %>%
    mutate(mean = binom$mean, l=binom$lower, u=binom$upper) %>%
    select(gs, phase, x, n, mean, l, u) -> rs_long
  
  denoms_structure = tibble(gs=c('no','yes'))
  
  long %>%
    mutate(gs = as.character(gs)) %>%
    filter(phase != 'Preclinical') %>%
    group_by(gs) %>%
    summarize(.groups='keep',
              denom = length(unique(ti_uid))) %>%
    ungroup() %>%
    right_join(denoms_structure, by='gs') %>%
    mutate(denom = replace_na(denom, 0)) -> denoms
  
  rs_long %>%
    filter(phase != 'Preclinical' & phase != 'I-Launch') %>%
    group_by(gs) %>%
    summarize(.groups='keep',
              x = x[phase=='III'],
              m = prod(mean),
              l = prod(mean) - qnorm(1 - (alpha)/2) * sqrt(prod(mean * (1 - mean) / n + mean^2) - prod(mean)^2), 
              u = prod(mean) + qnorm(1 - (alpha)/2) * sqrt(prod(mean * (1 - mean) / n + mean^2) - prod(mean)^2)) %>%
    rename(mean=m) %>%
    ungroup() %>%
    left_join(denoms, by='gs') %>%
    mutate(n = denom) %>%
    mutate(phase='I-Launch') %>%
    select(gs, phase, x, n, mean, l, u) -> ilrows
  
  
  rs_long %>%
    bind_rows(ilrows) %>%
    pivot_wider(id_cols=phase, names_from = gs, names_sep='_', values_from=c(x,n,mean,l,u)) %>%
    inner_join(select(phase_map, phase, phorder), by='phase') %>%
    arrange(phorder) %>%
    select(-phorder) %>%
    ungroup() %>%
    mutate(
      binratio = suppressWarnings(binom_ratio(x_yes, n_yes, x_no, n_no, mean_yes, mean_no, alpha = 0.05))) %>%
    mutate(
      rs_mean = binratio$rs_mean,
      rs_l = binratio$rs_l,
      rs_u = binratio$rs_u,
      fraction = glue("{x_yes}/{n_yes}") 
    ) %>%
    select(phase, x_yes, n_yes, x_no, n_no, mean_yes, l_yes, u_yes, mean_no, l_no, u_no, rs_mean, rs_l, rs_u, fraction) -> rr
  
  return(rr)
}
