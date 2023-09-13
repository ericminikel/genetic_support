

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
