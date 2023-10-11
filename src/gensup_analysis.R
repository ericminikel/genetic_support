overall_start_time = Sys.time()
cat(file=stderr(), 'Loading dependencies...')

options(stringsAsFactors=F)
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(binom))
suppressMessages(library(glue))
suppressMessages(library(lawstat))
suppressMessages(library(weights))
suppressMessages(library(epitools))
suppressMessages(library(DescTools))
suppressMessages(library(openxlsx))
suppressMessages(library(optparse))
suppressMessages(library(MASS)); summarize=dplyr::summarize; select=dplyr::select;
if(interactive()) {
  setwd('~/d/sci/src/genetic_support')
}

option_list = list(
  make_option(c("-o", "--oto"), action="store_true", default=FALSE, 
              help="limit to drugs with one target only (oto mode) [default %default]")
); 
opt = parse_args(OptionParser(option_list=option_list))

output_path = ifelse(opt$oto, 'oto/', 'display_items/')

##############
# OUTPUT FILE
##############

cat(file=stderr(), 'done.\nCreating output streams...'); flush.console()

text_stats_path = paste0(output_path,'stats_for_text.txt')
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T

supplement_path = paste0(output_path, 'supplement.xlsx')
supplement = createWorkbook()
supplement_directory = tibble(name=character(0), title=character(0))
write_supp_table = function(tbl, title='') {
  # write Excel sheet for supplement
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(tbl,paste0(output_path,'table_',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}



##############
# DATA INPUTS
##############

cat(file=stderr(), 'done.\nReading in data...')

merge2 = read_tsv('data/merge2.tsv.gz', col_types=cols())
pp = read_tsv('data/pp.tsv', col_types=cols())
drug_phase_summary = read_tsv('data/drug_phase_summary.tsv', col_types=cols())
assoc = read_tsv('data/assoc.tsv.gz', col_types=cols())
indic = read_tsv("data/indic.tsv", col_types=cols())
indic_topl_match = read_tsv('data/indic_topl_match.tsv', col_types=cols())
universe = read_tsv('data/universe.tsv', col_types=cols())
meta_hcat = read_tsv('data/meta_hcat.tsv', col_types=cols())
meta_acat = read_tsv('data/meta_acat.tsv', col_types=cols())
meta_ccat = read_tsv('data/meta_ccat.tsv', col_types=cols())
mesh_best_names = read_tsv('data/mesh_best_names.tsv.gz', col_types=cols())
sim = read_tsv('data/sim.tsv.gz', col_types=cols())

if (opt$oto) {
  pp$hcat = pp$oto_hcat
  pp$acat = pp$oto_acat
  pp$ccat = pp$oto_ccat
  merge2$hcat = pp$oto_hcat[match(merge2$ti_uid, pp$ti_uid)]
  merge2$acat = pp$oto_acat[match(merge2$ti_uid, pp$ti_uid)]
  merge2$ccat = pp$oto_ccat[match(merge2$ti_uid, pp$ti_uid)]
  pp = pp[!is.na(pp$ccat),]
  merge2 = merge2[!is.na(merge2$ccat),]
}

# constants

active_clinical = tibble(cat=c('Phase I','Phase II','Phase III'))

####
# Functions
####

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x < 0, '-', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}
lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

abs_or = function(odds_ratio) {
  abs_odds_ratio = odds_ratio
  flip_indices = odds_ratio < 1 & !is.na(odds_ratio)
  abs_odds_ratio[flip_indices] = 1/odds_ratio[flip_indices]
  return (abs_odds_ratio)
}

pipeline_best = function(merged_table,
                         basis='ti',
                         phase='historical',
                         require_insight=TRUE,
                         share_mode='L2G', # other option is V2G
                         min_share=0.5,
                         max_share=1,
                         worst_rank=Inf, # set to Inf if you want to include all
                         min_h4 = 0.9,
                         include_missing=FALSE,
                         associations=c('OMIM','GWAS'),
                         otg_subcat=c(''),
                         mendelian_mechanism='',
                         min_year=2005,
                         max_year=2022,
                         firstyear=F,
                         minusomim=F,
                         minusothersubcat=F,
                         mingenecount=0,
                         maxgenecount=Inf,
                         mapping_basis='all',
                         min_beta=0,
                         max_beta=Inf,
                         min_or=1,
                         max_or=Inf,
                         min_maf=0,
                         max_maf=1,
                         threshold=0.8,
                         network_list=NA,
                         verbose=T) {
  
  start_time = Sys.time()
  
  mtable = merged_table
  if (verbose) {
    cat(file=stderr(),'Starting row count: ',nrow(mtable),'\n')
    flush.console()
  }
  
  # add & select unique ID
  if (basis %in% c('di_mesh','drug-indication')) {
    mtable$di_uid = paste0(mtable$drugid,'-',mtable$indication_mesh_id)
    mtable$uid = mtable$di_uid
  } else if (basis %in% c('ti','target-indication')) {
    mtable$ti_uid = paste0(mtable$gene,'-',mtable$indication_mesh_id)
    mtable$uid = mtable$ti_uid
  } else if (basis=='drug') {
    mtable$uid = mtable$drugid
  }
  
  # assign highest level of advancement depending on phase specified
  if (phase == 'active') {
    meta = meta_acat
    mtable$cat = mtable$acat
  } else if (phase == 'historical') {
    meta = meta_hcat
    mtable$cat = mtable$hcat
  } else if (phase == 'combined') {
    meta = meta_ccat
    mtable$cat = mtable$ccat
  }
  
  mtable$catnum = meta$num[match(mtable$cat, meta$cat)]
  # remove "Other"
  mtable = mtable[mtable$cat != 'Other' ,]
  # map L2G share
  mtable$assoc_share = mtable$l2g_share
  mtable$assoc_rank = mtable$l2g_rank
  
  if (verbose) {
    cat(file=stderr(),'Selecting user-specified filters...')
    flush.console()
  }
  
  # by default, require non-missing target & indication
  if (!include_missing) {
    mtable = mtable[mtable$gene != '' & mtable$indication_mesh_id != '',]
  }
  # genetic insight requirement
  if (require_insight) {
    mtable = mtable[mtable$indication_mesh_id %in% indic$indication_mesh_id[indic$genetic_insight != 'none'],]
  }
  
  # remove omim-supported associations if desired. only works in T-I mode
  # note that order of operations is important - this must come before associations filter
  if (minusomim) {
    # look for first year in which a target-*indication* pair was genetically supported
    mtable %>%
      filter(comb_norm >= threshold) %>%
      filter(assoc_source %in% c('OMIM')) %>%
      select(ti_uid) -> omim_supported_ti
    # retain the null rows (i.e. no association) or those where hte T-I is not in OMIM
    # what gets removed? e.g. OTG associations that were already established by OMIM
    mtable %>%
      filter(is.na(assoc_source) | !(ti_uid %in% omim_supported_ti$ti_uid)) -> mtable
  }
  
  # user-specified sources of genetic associations
  # allow user to specify either grouping terms like "GWAS", or specific sources
  source_map = tibble(source=c("OTG", "PICCOLO", "Genebass", "OMIM", "intOGen"),
                      source_name=c('GWAS','GWAS','GWAS','OMIM','Somatic'))
  if (!(identical(associations, c('OMIM','GWAS','Somatic'))) ) {
    mtable %>%
      left_join(source_map, by=c('assoc_source'='source')) %>%
      filter(is.na(source_name) | source_name %in% associations | assoc_source %in% associations) -> mtable
    
  }
  
  # further filter of subtype of OTG association
  if (otg_subcat != '') {
    # first subset to just OTG
    mtable %>%
      filter(is.na(assoc_source) | assoc_source %in% 'OTG') -> mtable
    # now pick the types
    mtable %>%
      mutate(gwas_source = case_when(grepl('GCST', original_link) ~ 'GWAS Catalog',
                                     grepl('FINNGEN', original_link) ~ 'FinnGen',
                                     grepl('NEALE',original_link) ~ 'Neale UKBB',
                                     TRUE ~ 'Other')) %>%
      filter(is.na(assoc_source) | gwas_source %in% otg_subcat) -> mtable
  }
  
  
  # apply user-specified filter of OMIM disease mechanism
  if (mendelian_mechanism != '') {
    mtable %>%
      filter(!(mtable$assoc_source %in% 'OMIM') | is.na(mtable$assoc_info) | grepl(mendelian_mechanism,mtable$assoc_info)) -> mtable
  }
  
  # apply user-specified OTG gene mapping share & rank minimum/maximum
  if (share_mode == 'V2G') {
    mtable$assoc_share = mtable$v2g_share
    mtable$assoc_rank = mtable$v2g_rank
    assoc$assoc_share = assoc$v2g_share # needed in assoc table too for genecount section below
  } else if (share_mode == 'L2G') {
    mtable$assoc_share = mtable$l2g_share
    mtable$assoc_rank = mtable$l2g_rank
    assoc$assoc_share = assoc$l2g_share
  }
  
  # worst rank
  if (worst_rank < Inf) {
    mtable %>%
      filter(!(assoc_source %in% 'OTG') | mtable$assoc_rank <= worst_rank) -> mtable
  }
  
  # note that among OTG associations, throw out any with NA share, as these would be zeroes (does not occur in Dec 2021 dataset anyway)
  # and note that with L2G a significant number of associations actually have 100% share, so only delete > max_share and not >= max_share
  mtable %>%
    filter(!(assoc_source %in% 'OTG') | (!is.na(assoc_share) & assoc_share >= min_share & assoc_share <= max_share)) -> mtable
  
  # apply user-specified H4 minimum / maximum
  if (min_h4 > .9) {
    mtable %>%
      filter(!(assoc_source %in% 'PICCOLO') | (!is.na(mtable$pic_h4) & mtable$pic_h4 >= min_h4)) -> mtable
  }
  
  # apply user-specified genecount minimum/maximum
  if (mingenecount > 0 | maxgenecount < Inf) {
    assoc %>%
      filter(source %in% associations) %>%
      filter(source!='OTG' | (assoc_share >= min_share & assoc_share <= max_share)) %>%
      group_by(mesh_id) %>%
      summarize(.groups='keep', n_genes=length(unique(gene))) -> gene_counts
    mtable$gene_count = gene_counts$n_genes[match(mtable$assoc_mesh_id, gene_counts$mesh_id)]
    mtable %>%
      filter(is.na(gene_count) | gene_count >= mingenecount & gene_count <= maxgenecount) -> mtable
  }
  
  
  # apply "first year" criteria if applicable
  if (firstyear) {
    # look for first year in which a target-*indication* pair was genetically supported
    # only works in T-I mode
    mtable %>%
      filter(comb_norm >= threshold) %>%
      group_by(ti_uid) %>%
      summarize(.groups='keep', min_assoc_year=min(assoc_year)) -> ti_first_sup
    mtable$min_assoc_year = ti_first_sup$min_assoc_year[match(mtable$ti_uid, ti_first_sup$ti_uid)]
    # keep entries with no assoc year (the null rows), or where assoc year = the min assoc year, i.e. this is
    # the first report of this genetic association (or tied for first)
    mtable %>%
      filter(is.na(assoc_year) | assoc_year == min_assoc_year) -> mtable
  }
  
  # avoid -Inf values in comparisons by hard coding in case of all missing values:
  if (sum(!is.na(mtable$assoc_year)) == 0) {
    mtable_intrinsic_max_year = 2021
    mtable_intrinsic_min_year = 2000
  } else {
    mtable_intrinsic_max_year = max(mtable$assoc_year, na.rm=T)
    mtable_intrinsic_min_year = min(mtable$assoc_year, na.rm=T)
  }
  
  # apply user-specified filter of association years - for OTG only
  if (min_year > mtable_intrinsic_min_year | max_year < mtable_intrinsic_max_year) {
    mtable %>%
      filter(is.na(assoc_source) | assoc_source != 'OTG' | assoc_year >= min_year & assoc_year <= max_year) -> mtable
  }
  
  # join back in beta
  mtable$abs_beta = abs(assoc$beta[match(mtable$arow, assoc$arow)])
  if (min_beta > 0 | max_beta < Inf) {
    stopifnot(associations=='OTG') # only supported for OTG-only mode
    mtable %>%
      filter(is.na(assoc_source) | (!is.na(mtable$abs_beta) & mtable$abs_beta >= min_beta & mtable$abs_beta <= max_beta)) -> mtable
  }
  
  # same as beta but for OR
  mtable$abs_or = abs_or(assoc$odds_ratio[match(mtable$arow, assoc$arow)])  
  if (min_or > 1 | max_or < Inf) {
    stopifnot(associations=='OTG') # only supported for OTG-only mode
    # leave null rows (with no association source) but delete all those mapped to OTG that do not have or, or have or outside range
    mtable %>%
      filter(is.na(assoc_source) | (!is.na(mtable$abs_or) & mtable$abs_or >= min_or & mtable$abs_or <= max_or)) -> mtable
  }
  
  
  # lead SNP maf
  if (min_maf > 0 | max_maf < 1) {
    mtable$lead_maf = pmin(mtable$af_gnomad_nfe, 1-mtable$af_gnomad_nfe)
    mtable$lead_maf[!is.na(mtable$lead_maf) & mtable$lead_maf < 0] = NA
    # lead_maf >= min_maf & lead_maf < max_maf
    # >= and < gets you "[, )" logic
    # also remove those that are GWAS where lead_maf is NA - likely in non-European populations so af_gnomad_nfe is not relevant
    mtable %>% 
      filter(is.na(assoc_source) | !(assoc_source %in% c('OTG','PICCOLO','Genebass')) | (!is.na(lead_maf) & (lead_maf >= min_maf & lead_maf < max_maf))) -> mtable 
  }
  
  
  
  
  if (verbose) {
    cat(file=stderr(),'Selecting highest phase reached and best genetic similarity...')
    flush.console()
  }
  
  suppressWarnings(mtable %>% group_by(uid) %>% summarize(maxsim = max(comb_norm, na.rm=T), maxcat=max(catnum, na.rm=T)) -> step1)
  
  if (verbose) {
    cat(file=stderr(),nrow(step1),'rows remain.\n')
    flush.console()
  }
  
  if (verbose) {
    cat(file=stderr(),'Joining back in program details...')
    flush.console()
  }
  # add a filter first - only slightly reduces row count
  mtable %>%
    filter(uid %in% step1$uid & comb_norm %in% unique(step1$maxsim) & catnum %in% unique(step1$maxcat)) -> mtable
  # use tidy to left join
  step1 %>%
    left_join(mtable, by = c("uid" = "uid", "maxsim" = "comb_norm", "maxcat" = "catnum")) %>%
    rename(similarity=maxsim, catnum=maxcat) -> step2
  
  if (verbose) {
    cat(file=stderr(),nrow(step2),'rows found after join.\n')
    flush.console()
  }
  if (verbose) {
    cat(file=stderr(),'Removing duplicates resulting from ties...')
    flush.console()
  }
  
  # prioritize rows with known outcome at most advanced phase, then de-dup
  step2 %>%
    mutate(highest_phase_with_known_outcome = case_when(!is.na(succ_3_a) ~ 3,
                                                        !is.na(succ_2_3) ~ 2,
                                                        !is.na(succ_1_2) ~ 1,
                                                        !is.na(succ_p_1) ~ 0)) %>%
    arrange(uid, desc(highest_phase_with_known_outcome)) %>%
    group_by(uid) %>%
    slice(1) %>%
    ungroup() -> step2 # row count should drop back to ~ that of step1

  # this is how we used to accomplish the above, with full dataset - drug ids and years
  # if (basis=='drug') {
  #   # if a drug has both active and historical programs, prioritize historical since it has a known outcome
  #   step2 = step2[with(step2, order(uid, -year_use_valid, -hcatnum)),]
  # } else {
  #   # simply prioritize the programs with valid years for current phase
  #   step2 = step2[with(step2, order(uid, -year_use_valid)),] 
  # }
  
  
  if (verbose) {
    cat(file=stderr(),nrow(step2),'rows remain.\n')
    flush.console()
  }
  
  # annotate in additional info:
  step2$areas = indic$areas[match(step2$indication_mesh_id, indic$indication_mesh_id)]
  step2$genetic_insight = replace_na(indic$genetic_insight[match(step2$indication_mesh_id, indic$indication_mesh_id)],'none')
  step2$target_status = ''
  step2$target_status[step2$similarity >= threshold] = 'genetically supported target'
  step2$target_status[step2$similarity <  threshold] = 'unsupported target'
  step2$target_status[step2$genetic_insight == 'none'] = 'indication lacks genetic insight'
  step2$target_status[is.na(step2$gene) | step2$gene == ''] = 'no target annotated'
  step2$target_status[is.na(step2$indication_mesh_id) | step2$indication_mesh_id == ''] = 'no indication annotated'
  

  if (verbose) {
    cat(file=stderr(),paste0('Using sim threshold ',threshold,', "genetically supported target" rows: ',sum(step2$target_status=='genetically supported target'),'....\n'))
    time_elapsed = (Sys.time() - start_time)
    cat(file=stderr(),'pipeline_best completed in',round(time_elapsed,1),units(time_elapsed),'.\n')
    flush.console()
  }
  
  return (step2)
}


advancement_forest = function(best_table, phase='historical') {
  if (phase == 'active') {
    meta = meta_acat
  } else if (phase == 'historical') {
    meta = meta_hcat
  } else if (phase == 'combined') {
    meta = meta_ccat
  }
  meta %>%
    left_join(best_table, by=c('num'='catnum', 'cat'='cat')) %>%
    filter(cat != 'Other') %>%
    filter(!(target_status %in% c('indication lacks genetic insight','no indication annotated','no target annotated'))) %>%
    group_by(catnum=num, cat) %>%
    summarize(.groups='keep',
              n_total = n(),
              n_gensup = sum(target_status=='genetically supported target', na.rm=T)) -> forest_data
  bconf_obj = binom.confint(x=forest_data$n_gensup, n=forest_data$n_total, method='wilson', conf.level = .95)
  forest_data$proportion = bconf_obj$mean
  forest_data$l95 = bconf_obj$lower
  forest_data$u95 = bconf_obj$upper
  forest_data$y = max(forest_data$catnum) - forest_data$catnum + 1
  colnames(forest_data) = c('num','label','denominator','numerator','mean','l95','u95','y')
  return (forest_data)
}

advancement_rr = function(best, alpha = 0.05, threshold = NA) {
  
  phase_map = tibble(phase=factor(c('Preclinical','I','II','III','I-Launch'),ordered=T,levels=c('Preclinical','I','II','III','I-Launch')), 
                     phorder = 0:4,
                     varname=c('succ_p_1','succ_1_2','succ_2_3','succ_3_a','succ_1_a'))

  if (is.na(threshold)) {
    best$gensup = best$target_status=='genetically supported target'
  } else {
    best$gensup = !is.na(best$similarity) & best$similarity >= threshold
  } 
  
  best %>%
    filter(genetic_insight != 'none') %>%
    select(ti_uid, gensup, succ_p_1:succ_3_a) %>%
    pivot_longer(succ_p_1:succ_3_a) %>%
    inner_join(phase_map, by=c('name'='varname')) %>%
    filter(!is.na(value)) %>%
    rename(success=value) %>%
    mutate(gs = ifelse(gensup,'yes','no')) %>%
    select(ti_uid, gs, phase, success) -> long
  
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
  
  rr_long_structure = crossing(gs=c('yes','no'),phase=c('Preclinical','I','II','III'))
  
  # Wilson CI for single phases
  long %>%
    mutate(gs = as.character(gs)) %>%
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
  
  # Wald CI for product of P(S) across I-Launch
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
  
  return (rr)
}

always_katz = T
binom_ratio_atomic = function(x_yes, n_yes, x_no, n_no, mean_yes=NA, mean_no=NA, alpha =  0.05) {
  
  # TO DO: think through whether we really need the Katz option, esp since for 
  # I-Launch it will be based on the combined denominator (unique T-I in any phase)
  
  # Wald by default
  if (!always_katz & x_yes >= 3 & x_no >=3) {
    mean = mean_yes / mean_no
    #mean = (x_yes/n_yes) / (x_no/n_no)
    lower = mean - qnorm(1 - alpha/2) * sqrt(1/(n_yes + n_no) * mean^2 * (1/mean_yes + 1/mean_no))
    upper = mean + qnorm(1 - alpha/2) * sqrt(1/(n_yes + n_no) * mean^2 * (1/mean_yes + 1/mean_no))
  } else { # Katz for small N
    # for the I-Launch row, mean_yes or mean_no may be NA because a phase had no data
    # important to leave that NA - don't use the filled-in total denominator
    if (is.na(mean_yes) | is.na(mean_no)) {
      mean = lower = upper = as.numeric(NA)
    } else {
      binom_obj = BinomRatioCI(x_yes, n_yes, x_no, n_no, conf=1-alpha)
      mean = binom_obj[1,'est']
      lower = binom_obj[1,'lwr.ci']
      upper = binom_obj[1,'upr.ci']
    }
  }
  return (cbind(mean=mean,lower=lower,upper=upper))
}

binom_ratio = function(x_yes, n_yes, x_no, n_no, mean_yes, mean_no, alpha =  0.05) {
  setNames(as_tibble(t(mapply(binom_ratio_atomic, x_yes, n_yes, x_no, n_no, mean_yes, mean_no, alpha))), c('rs_mean','rs_l','rs_u'))
}


plot_forest = function(forestdf, xlims=c(0,1), xstyle='percent', mar=c(3,8,3,8), xlab='', title='', col='#000000', showvals=F, right_text=NA, xlab_line=1.6) {
  ylims = range(forestdf$y) + c(-0.5, 0.5)
  par(mar=mar)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  if (xstyle == 'percent') {
    axis(side=1, at=0:100/100, labels=NA, tck=-0.025)
    if (max(xlims) > .5) {
      axis(side=1, at=0:10/10, labels=NA, tck=-0.05)
      axis(side=1, at=0:2/2, labels=percent(0:2/2), lwd=0, line=-0.5)
    } else {
      axis(side=1, at=0:20/20, labels=NA, tck=-0.05)
      axis(side=1, at=0:20/20, labels=percent(0:20/20,digits=0), lwd=0, line=-0.5)
    }
  } else if (xstyle == 'ratio') {
    axis(side=1, at=seq(0,max(xlims),0.1), labels=NA, tck=-0.025)
    axis(side=1, at=seq(0,max(xlims),1), labels=NA, tck=-0.05)
    axis(side=1, at=seq(0,max(xlims),1), lwd=0, line=-0.5)
    abline(v=1, lwd=0.25, lty=3)
  }
  mtext(side=1, line=xlab_line, cex=0.75, text=xlab)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  mtext(side=2, at=forestdf$y, text=forestdf$label, cex=0.75, line=0.5, las=2, col=col)
  mtext(side=4, at=forestdf$y, text=paste0(formatC(forestdf$numerator,big.mark=','),'/',formatC(forestdf$denominator,big.mark=',')), cex=0.75, las=2, line=0.25)
  par(xpd=T)
  if (is.na(right_text)) {
    if (xstyle=='ratio') {
      right_text = 'Approved/\nSupported'
    } else if (xstyle == 'percent') {
      right_text = 'Supported/\nTotal'
    }
  }
  mtext(side=4, las=2, at=max(forestdf$y)+1, font=2, text=right_text, cex=0.7, padj=0)
  par(xpd=F)
  points(forestdf$mean, forestdf$y, pch=19, col=col) # means
  segments(x0=forestdf$l95, x1=pmin(forestdf$u95,max(xlims)), y0=forestdf$y, lwd=2, col=col) # 95%CIs
  mtext(side=3, line=0, text=title, col=col, font=1, cex=0.7)
  if (showvals) {
    text(x=forestdf$u95, y=forestdf$y, pos=4, labels=formatC(forestdf$mean, format='f', digits=1), font=3, cex=.75, col=col)
  }
}


subset_by_area = function(best_table, topl, filter='only', orphan='any') {
  btable = best_table
  if (orphan=='any') {
    btable = btable
  } else if (orphan=='yes') {
    btable = btable[btable$orphan==1,]
  } else if (orphan=='non') {
    btable = btable[btable$orphan==0,]
  }
  if (topl == 'all' | topl == 'ALL') {
    btable = btable
  } else if (filter == 'only') {
    btable = btable[btable$indication_mesh_id %in% indic_topl_match$indication_mesh_id[indic_topl_match$topl==topl],]
  } else if (filter == 'non') {
    btable = btable[!(btable$indication_mesh_id %in% indic_topl_match$indication_mesh_id[indic_topl_match$topl==topl]),]
  }
  
  return (btable)
}


########
# Staging
########

cat(file=stderr(), 'done.\nGenerating pipeline tables: hist_ti...')
hist_ti = pipeline_best(merge2, phase='historical', basis='ti', verbose = F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_all...')
hist_ti_all = pipeline_best(merge2, phase='historical', basis='ti', require_insight=F, include_missing=F, verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: active_ti...')
active_ti = pipeline_best(merge2, phase='active', basis='ti', verbose=F)

cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_everything...')
hist_ti_everything = pipeline_best(merge2, phase='historical', basis='ti', associations = c('OMIM','GWAS','intOGen'), verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_germline...')
hist_ti_germline = hist_ti
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_omim...')
hist_ti_omim = pipeline_best(merge2, phase='historical', basis='ti', associations=c('OMIM'), verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_genebass...')
hist_ti_genebass = pipeline_best(merge2, phase='historical', basis='ti', associations=c('Genebass'), verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_otg...')
hist_ti_otg = pipeline_best(merge2, phase='historical', basis='ti', associations=c('OTG'), verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_piccolo...')
hist_ti_pic  = pipeline_best(merge2, phase='historical', basis='ti', associations=c('PICCOLO'), verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_gwas...')
hist_ti_gwas  = pipeline_best(merge2, phase='historical', basis='ti', associations=c('PICCOLO','OTG','Genebass'), verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_soma...')
hist_ti_soma = pipeline_best(merge2, phase='historical', basis='ti', associations=c('Somatic'), verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_intogen...')
hist_ti_intogen = pipeline_best(merge2, phase='historical', basis='ti', associations=c('intOGen'), verbose=F)

cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_gwascat...')
hist_ti_gwascat  = pipeline_best(merge2, phase='historical', basis='ti', otg_subcat='GWAS Catalog', verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_ukbb...')
hist_ti_ukbb  = pipeline_best(merge2, phase='historical', basis='ti', otg_subcat='Neale UKBB', verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: hist_ti_finngen...')
hist_ti_finngen  = pipeline_best(merge2, phase='historical', basis='ti', otg_subcat='FinnGen', verbose=F)

cat(file=stderr(), '\rGenerating pipeline tables: combined_ti...')
combined_ti = pipeline_best(merge2,  phase='combined', basis='ti', verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: combined_ti_omim...')
combined_ti_omim = pipeline_best(merge2,  phase='combined', basis='ti', associations=c('OMIM'), verbose=F)
cat(file=stderr(), '\rGenerating pipeline tables: combined_ti_gwas...')
combined_ti_gwas = pipeline_best(merge2,  phase='combined', basis='ti', associations=c('OTG','PICCOLO','Genebass'), verbose=F)

cat(file=stderr(), '\rGenerating pipeline tables... done.                    \n')

# combined_ti %>%
#   select(ti_uid, gene, indication_mesh_id, indication_mesh_term, ccatnum, ccat, orphan, 
#          assoc_mesh_id, assoc_mesh_term, assoc_source, assoc_info, original_trait, original_link,
#          assoc_year, pic_qtl_pval, pic_h4)


#########
# Figure S1
#########

# Things I can't generate based on my limited dataset:
# Unique drugs - total
# Unique drugs - monotherapy, phase assigned, human target

merge2$tia = paste(merge2$gene, merge2$indication_mesh_id, merge2$assoc_mesh_id, sep='-')
merge2$tia[is.na(merge2$gene) | is.na(merge2$indication_mesh_id) | is.na(merge2$assoc_mesh_id)] = NA

write(paste('Unique indications: ',nrow(indic),'\n',sep=''),text_stats_path,append=T)
write(paste('Unique indications - of which genetic insight: ',sum(indic$genetic_insight!='none'),'\n',sep=''),text_stats_path,append=T)
write(paste('Unique targets: ',length(unique(pp$gene)),'\n',sep=''),text_stats_path,append=T)
write(paste('Unique T-I: ',nrow(pp),'\n',sep=''),text_stats_path,append=T)
write(paste('Unique D-T-I: ',nrow(drug_phase_summary),'\n',sep=''),text_stats_path,append=T)

write(paste('Approved | Unique targets: ',length(unique(pp$gene[pp$ccat=='Launched'])),'\n',sep=''),text_stats_path,append=T)
write(paste('Approved | Unique indications: ',length(unique(pp$indication_mesh_id[pp$ccat=='Launched'])),'\n',sep=''),text_stats_path,append=T)
write(paste('Approved | Unique T-I: ',sum(pp$ccat=='Launched'),'\n',sep=''),text_stats_path,append=T)

write(paste('Historical | Unique targets: ',length(unique(pp$gene[pp$hcat!='Launched'])),'\n',sep=''),text_stats_path,append=T)
write(paste('Historical | Unique indications: ',length(unique(pp$indication_mesh_id[pp$hcat!='Launched'])),'\n',sep=''),text_stats_path,append=T)
write(paste('Historical | Unique T-I: ',sum(pp$hcat!='Launched'),'\n',sep=''),text_stats_path,append=T)

write(paste('Active | Unique targets: ',length(unique(pp$gene[pp$acat!='Launched' & !is.na(pp$acat)])),'\n',sep=''),text_stats_path,append=T)
write(paste('Active | Unique indications: ',length(unique(pp$indication_mesh_id[pp$acat!='Launched' & !is.na(pp$acat)])),'\n',sep=''),text_stats_path,append=T)
write(paste('Active | Unique T-I: ',sum(pp$acat!='Launched' & !is.na(pp$acat)),'\n',sep=''),text_stats_path,append=T)

write(paste('Merged sim ≥0.8 | Unique targets: ',length(unique(merge2$gene[merge2$comb_norm >= 0.8 & !is.na(merge2$comb_norm)])),'\n',sep=''),text_stats_path,append=T)
write(paste('Merged sim ≥0.8 | Unique indications: ',length(unique(merge2$indication_mesh_id[merge2$comb_norm >= 0.8 & !is.na(merge2$comb_norm)])),'\n',sep=''),text_stats_path,append=T)
write(paste('Merged sim ≥0.8 | Unique T-I: ',length(unique(merge2$ti_uid[merge2$comb_norm >= 0.8 & !is.na(merge2$comb_norm)])),'\n',sep=''),text_stats_path,append=T)
write(paste('Merged sim ≥0.8 | Unique T-I-A: ',length(unique(merge2$tia[merge2$comb_norm >= 0.8 & !is.na(merge2$comb_norm)])),'\n',sep=''),text_stats_path,append=T)

assoc$ta_uid = paste0(assoc$gene, '-', assoc$mesh_id)
write(paste('Assocs | Unique T-A: ',length(unique(assoc$ta_uid[(assoc$l2g_share >= 0.5 | assoc$source != 'OTG')])),'\n',sep=''),text_stats_path,append=T)
write(paste('Assocs | OTG | Unique T-A: ',length(unique(assoc$ta_uid[(assoc$l2g_share >= 0.5 & assoc$source == 'OTG')])),'\n',sep=''),text_stats_path,append=T)
write(paste('Assocs | OMIM | Unique T-A: ',length(unique(assoc$ta_uid[(assoc$source == 'OMIM')])),'\n',sep=''),text_stats_path,append=T)
write(paste('Assocs | PICCOLO | Unique T-A: ',length(unique(assoc$ta_uid[(assoc$source == 'PICCOLO')])),'\n',sep=''),text_stats_path,append=T)
write(paste('Assocs | Genebass | Unique T-A: ',length(unique(assoc$ta_uid[(assoc$source == 'Genebass')])),'\n',sep=''),text_stats_path,append=T)
write(paste('Assocs | intOGen | Unique T-A: ',length(unique(assoc$ta_uid[(assoc$source == 'intOGen')])),'\n',sep=''),text_stats_path,append=T)

assoc %>%
  filter(!is.na(source)) %>%
  filter(!is.na(mesh_id)) %>%
  distinct(mesh_id) %>%
  mutate(in_sim = mesh_id %in% sim$meshcode_a) -> assoc_mesh_uq
pp %>%
  distinct(indication_mesh_id) %>%
  filter(!is.na(indication_mesh_id)) %>%
  mutate(in_sim = indication_mesh_id %in% sim$meshcode_a) -> pp_mesh_uq

write(paste0('Pharmaprojects | MeSH missing from sim matrix: ',sum(!(pp_mesh_uq$indication_mesh_id %in% sim$meshcode_a)),'/',nrow(pp_mesh_uq),'\n',sep=''),text_stats_path,append=T)
write(paste0('Pharmaprojects | MeSH present in sim matrix: ',sum((pp_mesh_uq$indication_mesh_id %in% sim$meshcode_a)),'/',nrow(pp_mesh_uq),'\n',sep=''),text_stats_path,append=T)
write(paste0('Pharmaprojects | proportion of unique MeSH in sim matrix: ',percent(mean(pp_mesh_uq$indication_mesh_id %in% sim$meshcode_a), digits=3),'\n',sep=''),text_stats_path,append=T)
write(paste0('Pharmaprojects | proportion assoc rows in sim matrix: ',percent(mean(pp$indication_mesh_id[!is.na(pp$indication_mesh_id)] %in% sim$meshcode_a), digits=3),'\n',sep=''),text_stats_path,append=T)

write(paste0('Assocs | MeSH missing from sim matrix: ',sum(!(assoc_mesh_uq$mesh_id %in% sim$meshcode_a)),'/',nrow(assoc_mesh_uq),'\n',sep=''),text_stats_path,append=T)
write(paste0('Assocs | MeSH present in sim matrix: ',sum((assoc_mesh_uq$mesh_id %in% sim$meshcode_a)),'/',nrow(assoc_mesh_uq),'\n',sep=''),text_stats_path,append=T)
write(paste0('Assocs | proportion of unique MeSH in sim matrix: ',percent(mean(assoc_mesh_uq$mesh_id %in% sim$meshcode_a), digits=3),'\n',sep=''),text_stats_path,append=T)
write(paste0('Assocs | proportion assoc rows in sim matrix: ',percent(mean(assoc$mesh_id[!is.na(assoc$source)] %in% sim$meshcode_a), digits=3),'\n',sep=''),text_stats_path,append=T)


#########
# Figure 1
#########

#### 1A
hist_ti_forest = advancement_forest(hist_ti,phase='historical')
active_ti_forest = advancement_forest(active_ti,phase='active')
combined_ti_forest = advancement_forest(combined_ti,phase='combined')

#### 1B
assoc_source_rr_forest = tibble(label=c('All germline','OMIM','All GWAS','All OTG','GWAS Catalog','Neale UKBB','FinnGen','PICCOLO','Genebass'),
                                pipeline_obj=c('hist_ti_germline','hist_ti_omim','hist_ti_gwas','hist_ti_otg','hist_ti_gwascat','hist_ti_ukbb','hist_ti_finngen','hist_ti_pic','hist_ti_genebass')) %>%
  mutate(y=max(row_number()) - row_number() + 1)
for (i in 1:nrow(assoc_source_rr_forest)) {
  pipeline_obj = get(assoc_source_rr_forest$pipeline_obj[i])
  rr_obj = advancement_rr(pipeline_obj)
  assoc_source_rr_forest[i,c('mean','l95','u95')] = rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  assoc_source_rr_forest[i,c('numerator')]        = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  assoc_source_rr_forest[i,c('denominator')]      = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
}

#### 1C
cat(file=stderr(), '\rCalculating ROC L2G...')
roc_l2g = data.frame(share=seq(0.1, 0.9, 0.05))
for (i in 1:nrow(roc_l2g)) {
  
  pipeline_obj = pipeline_best(merge2, phase='historical', basis='ti', associations=c('OTG'), share_mode='L2G', min_share=roc_l2g$share[i], verbose=F)
  
  if (sum(pipeline_obj$target_status=='genetically supported target')==0) {
    next
  }
  
  rr_obj = advancement_rr(pipeline_obj)
  roc_l2g[i,c('mean','l95','u95')] = rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  roc_l2g[i,c('numerator')]        = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  roc_l2g[i,c('denominator')]      = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
  roc_l2g[i,c('gensup_all')]       = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
  roc_l2g[i,c('gensup_launched')]  = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  roc_l2g[i,c('nosup_all')]        = rr_obj[rr_obj$phase=='I-Launch',c('n_no')]
  roc_l2g[i,c('nosup_launched')]   = rr_obj[rr_obj$phase=='I-Launch',c('x_no')]
  
}

#### 1D

cat(file=stderr(), 'done.\nCalculating year of discovery statistics...')

year_tranches = tibble(tranche=1:5,
                       minyear=c(2007,2007,2011,2015,2019),
                       maxyear=c(2022,2010,2014,2018,2022)) %>%
  mutate(label = paste(minyear,maxyear,sep='-')) %>%
  mutate(label = ifelse(tranche==1,'All',label))
modes = tibble(abbr = c('a','b','c'),
               firstyear=c(F,T,T),
               minusomim=c(F,F,T))
crossing(year_tranches, modes) %>%
  add_column(mean=NA, l95=NA, u95=NA, numerator=NA, denominator=NA) %>% 
  arrange(tranche) -> year_rrs

for (i in 1:nrow(year_rrs)) {
  pb_obj = pipeline_best(merge2, 
                         phase='historical', 
                         basis='ti', 
                         associations='OTG',
                         otg_subcat='GWAS Catalog',
                         min_year=year_rrs$minyear[i],
                         max_year=year_rrs$maxyear[i],
                         firstyear=year_rrs$firstyear[i],
                         minusomim=year_rrs$minusomim[i], verbose=F)
  rr_obj = advancement_rr(pb_obj)
  year_rrs$mean[i] = rr_obj$rs_mean[rr_obj$phase=='I-Launch']
  year_rrs$l95[i] = rr_obj$rs_l[rr_obj$phase=='I-Launch']
  year_rrs$u95[i] = rr_obj$rs_u[rr_obj$phase=='I-Launch']
  year_rrs$numerator[i] = rr_obj$x_yes[rr_obj$phase=='I-Launch']
  year_rrs$denominator[i] = rr_obj$n_yes[rr_obj$phase=='I-Launch']
}

year_rrs$y = rep(5:1,each=3)
year_rrs$frac = paste0(year_rrs$numerator,'/',year_rrs$denominator)

yearfirst_all_pbest = pipeline_best(merge2, 
                                    phase='historical', 
                                    basis='ti', 
                                    associations='OTG',
                                    otg_subcat='GWAS Catalog',
                                    min_year=2007,
                                    max_year=2022,
                                    firstyear=T,
                                    minusomim=T, verbose=F)
yearfirst_logit_data = subset(yearfirst_all_pbest, target_status=='genetically supported target' & cat %in% c('Launched',active_clinical$cat))
yearfirst_logit_data$launched = yearfirst_logit_data$cat=='Launched'
yearfirst_logit = glm(launched ~ assoc_year, data=yearfirst_logit_data, family='binomial')
yearfirst_logit_beta = summary(yearfirst_logit)$coefficients['assoc_year','Estimate']
yearfirst_logit_p = summary(yearfirst_logit)$coefficients['assoc_year','Pr(>|z|)']

write(paste('Logit model launched ~ assoc_year: beta = ',formatC(yearfirst_logit_beta,digits=2,format='g'),', P = ',formatC(yearfirst_logit_p,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)

yearfirst_logit_data %>% 
  group_by(assoc_year) %>%
  summarize(.groups='keep', n_launched=sum(ccat=='Launched'), n_i_l=sum(ccat %in% c('Launched',active_clinical$cat))) -> yearfirst_bins

yearfirst_bins$p_launched = yearfirst_bins$n_launched / yearfirst_bins$n_i_l
yearfirst_bins[,c('p_launched','l95','u95')] = binom.confint(x=yearfirst_bins$n_launched,n=yearfirst_bins$n_i_l,method='wilson')[,c('mean','lower','upper')]
yearfirst_bins = yearfirst_bins[yearfirst_bins$n_i_l > 0,]

cat(file=stderr(), 'done.\nCalculating gene count statistics...')

gc_tranches = tibble(tranche=1:6,
                     mingenecount=c(1,1,2,10,100,1000),
                     maxgenecount=c(Inf,1,9,99,999,Inf),
                     label = c('All','1','2-9','10-99','100-999','1,000+'))

gc_modes = tibble(abbr = c('a'),
                  associations = c('OTG'),
                  subcats = c('GWAS Catalog'))

crossing(gc_tranches, gc_modes) %>%
  add_column(mean=NA, l95=NA, u95=NA, numerator=NA, denominator=NA) %>% 
  arrange(tranche) -> gc_rrs


for (i in 1:nrow(gc_rrs)) {
  pb_obj = pipeline_best(merge2, 
                         phase='historical', 
                         basis='ti', 
                         associations=gc_rrs$associations[i],
                         otg_subcat=gc_rrs$subcats[i],
                         mingenecount=gc_rrs$mingenecount[i],
                         maxgenecount=gc_rrs$maxgenecount[i], verbose=F)
  if (i == 1) {
    gc_logit_data = subset(pb_obj,  target_status=='genetically supported target' & !is.na(indication_mesh_id) & cat %in% c('Launched',active_clinical$cat))
  }
  rr_obj = advancement_rr(pb_obj)
  gc_rrs$mean[i] = rr_obj$rs_mean[rr_obj$phase=='I-Launch']
  gc_rrs$l95[i] = rr_obj$rs_l[rr_obj$phase=='I-Launch']
  gc_rrs$u95[i] = rr_obj$rs_u[rr_obj$phase=='I-Launch']
  gc_rrs$numerator[i] = rr_obj$x_yes[rr_obj$phase=='I-Launch']
  gc_rrs$denominator[i] = rr_obj$n_yes[rr_obj$phase=='I-Launch']
}

gc_logit_data$launched = gc_logit_data$cat=='Launched'
genecount_logit = glm(launched ~ gene_count, data=gc_logit_data, family='binomial')
genecount_logit_beta = summary(genecount_logit)$coefficients['gene_count','Estimate']
genecount_logit_p = summary(genecount_logit)$coefficients['gene_count','Pr(>|z|)']

write(paste('Logit model launched ~ gene_count: beta = ',formatC(genecount_logit_beta,digits=2,format='g'),', P = ',formatC(genecount_logit_p,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)


gc_logit_data$top3areas = grepl('hematologic|respiratory|metabolic',gc_logit_data$areas)
gc_logit_data$gc100plus = gc_logit_data$gene_count >= 100

top3_gc100_ctable = table(gc_logit_data[,c('top3areas','gc100plus')])
fisher_obj = fisher.test(top3_gc100_ctable)

write(paste('Enrichment of top 3 highest-RS therapy areas for gene_count >= 100: OR = ',formatC(fisher_obj$estimate,digits=2,format='g'),', P = ',formatC(fisher_obj$p.value,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)



cat(file=stderr(), 'done.\nCalculating beta statistics...')


ideal_beta_breaks = quantile(abs(assoc$beta[assoc$source=='OTG' & grepl('GCST',assoc$original_link)]), 0:4/4, na.rm=T)
ideal_beta_breaks_disp = round(ideal_beta_breaks,digits=3)
beta_rrs = tibble(y=5:1, 
                  label=c('All with beta',
                          paste0('beta 0 - ',ideal_beta_breaks_disp[2]),
                          paste0('beta ',ideal_beta_breaks_disp[2],' - ',ideal_beta_breaks_disp[3]),
                          paste0('beta ',ideal_beta_breaks_disp[3],' - ',ideal_beta_breaks_disp[4]),
                          paste0('beta ',ideal_beta_breaks_disp[4],'+')),
                  min_beta = c(0,     0, ideal_beta_breaks[2:4]),
                  max_beta = c(1e100, ideal_beta_breaks[2:4], 1e100))
for (i in 1:nrow(beta_rrs)) {
  pb_obj = pipeline_best(merge2, phase='historical', basis='target-indication', associations='OTG', otg_subcat='GWAS Catalog', min_beta=beta_rrs$min_beta[i], max_beta=beta_rrs$max_beta[i], verbose=F)
  if (i == 1) {
    beta_logit_data = subset(pb_obj,  target_status=='genetically supported target' & !is.na(indication_mesh_id) & cat %in% c('Launched',active_clinical$cat))
  }
  rr_obj = advancement_rr(pb_obj)
  beta_rrs[i,c('mean','l95','u95')] = rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  beta_rrs[i,c('numerator')]        = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  beta_rrs[i,c('denominator')]      = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
}

beta_logit_data$launched = beta_logit_data$cat=='Launched'
beta_logit = glm(launched ~ abs_beta, data=beta_logit_data, family='binomial')
beta_logit_beta = summary(beta_logit)$coefficients['abs_beta','Estimate']
beta_logit_p = summary(beta_logit)$coefficients['abs_beta','Pr(>|z|)']


write(paste('Logit model launched ~ abs_beta: beta = ',formatC(beta_logit_beta,digits=2,format='g'),', P = ',formatC(beta_logit_p,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)


cat(file=stderr(), 'done.\nCalculating OR statistics...')


ideal_or_breaks = quantile(abs_or(assoc$odds_ratio[assoc$source=='OTG' & grepl('GCST',assoc$original_link)]), 0:4/4, na.rm=T)
ideal_or_breaks_disp = round(ideal_or_breaks,digits=3)
or_rrs = tibble(y=5:1, 
                label=c('All with OR',
                        paste0('OR 1 - ',ideal_or_breaks_disp[2]),
                        paste0('OR ',ideal_or_breaks_disp[2],' - ',ideal_or_breaks_disp[3]),
                        paste0('OR ',ideal_or_breaks_disp[3],' - ',ideal_or_breaks_disp[4]),
                        paste0('OR ',ideal_or_breaks_disp[4],'+')),
                min_or = c(1,     1, ideal_or_breaks[2:4]),
                max_or = c(1e100, ideal_or_breaks[2:4], 1e100))
for (i in 1:nrow(or_rrs)) {
  pb_obj = pipeline_best(merge2, phase='historical', basis='target-indication', associations='OTG', otg_subcat='GWAS Catalog', min_or=or_rrs$min_or[i], max_or=or_rrs$max_or[i], verbose=F)
  if (i == 1) {
    or_logit_data = subset(pb_obj,  target_status=='genetically supported target' & !is.na(indication_mesh_id) & cat %in% c('Launched',active_clinical$cat))
  }
  rr_obj = advancement_rr(pb_obj)
  or_rrs[i,c('mean','l95','u95')] = rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  or_rrs[i,c('numerator')]        = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  or_rrs[i,c('denominator')]      = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
}

or_logit_data$launched = or_logit_data$cat=='Launched'
oddsratio_logit = glm(launched ~ abs_or, data=or_logit_data, family='binomial')
oddsratio_logit_beta = summary(oddsratio_logit)$coefficients['abs_or','Estimate']
oddsratio_logit_p = summary(oddsratio_logit)$coefficients['abs_or','Pr(>|z|)']

write(paste('Logit model launched ~ abs_or: beta = ',formatC(oddsratio_logit_beta,digits=2,format='g'),', P = ',formatC(oddsratio_logit_p,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)

cat(file=stderr(), 'done.\nCalculating MAF statistics...')

gwascat = suppressMessages(read_tsv('../digap/data/gwascat/gwas_catalog-ancestry_r2022-02-02.tsv') %>% clean_names())
gwascat$acc_link = paste0('https://genetics.opentargets.org/study/',gwascat$study_accession)
gwascat$nfe = gwascat$broad_ancestral_category=='European' & gwascat$country_of_origin != 'Finland'
merge2$in_nfe = gwascat$nfe[match(merge2$original_link, gwascat$acc_link)]
merge2$af_gnomad_nfe[!merge2$in_nfe | is.na(merge2$in_nfe)] = NA

assoc$lead_maf = ifelse(assoc$af_gnomad_nfe > 0.5, 1-assoc$af_gnomad_nfe, assoc$af_gnomad_nfe)
assoc$lead_maf[assoc$lead_maf < 0] = NA
maf_or_cor = cor.test(assoc$lead_maf, ifelse(assoc$odds_ratio > 1, assoc$odds_ratio, 1/assoc$odds_ratio))
maf_beta_cor = cor.test(assoc$lead_maf, abs(assoc$beta))
maf_rrs = tibble(y=5:1, 
                 label=c('All with MAF','MAF 1 - 3%','MAF 3 - 10%','MAF 10 - 30%','MAF 30 - 50%'),
                 min_maf = c(0,     0.01,   0.03,  0.10, .30),
                 max_maf = c(0.501, 0.03,   0.10,  0.30, .501))
for (i in 1:nrow(maf_rrs)) {
  pb_obj = pipeline_best(merge2, phase='historical', basis='target-indication', associations='OTG', otg_subcat='GWAS Catalog', firstyear=T, min_maf=maf_rrs$min_maf[i], max_maf=maf_rrs$max_maf[i], verbose=F)
  if (i == 1) {
    maf_logit_data = subset(pb_obj,  target_status=='genetically supported target' & !is.na(indication_mesh_id) & cat %in% c('Launched',active_clinical$cat))
  }
  rr_obj = advancement_rr(pb_obj)
  maf_rrs[i,c('mean','l95','u95')] = rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  maf_rrs[i,c('numerator')]        = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  maf_rrs[i,c('denominator')]      = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
}
maf_logit_data$launched = maf_logit_data$cat=='Launched'
maf_logit = glm(launched ~ lead_maf, data=maf_logit_data, family='binomial')
maf_logit_beta = summary(maf_logit)$coefficients['lead_maf','Estimate']
maf_logit_p = summary(maf_logit)$coefficients['lead_maf','Pr(>|z|)']


write(paste('Logit model launched ~ lead_maf: beta = ',formatC(maf_logit_beta,digits=2,format='g'),', P = ',formatC(maf_logit_p,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)

master_cols = c('label','mean','l95','u95','numerator','denominator')
rbind(year_rrs[year_rrs$abbr=='c',master_cols], 
      gc_rrs[,master_cols], 
      beta_rrs[,master_cols], 
      or_rrs[,master_cols],
      maf_rrs[,master_cols]) -> master_forest

cat(file=stderr(), 'done.\nCreating Figure 1...')

resx=300
png(paste0(output_path,'/figure-1.png'),width=6.5*resx,height=5.5*resx,res=resx)

#layout_matrix = matrix(c(1,2,3,4,4,4),nrow=3,byrow=F)
layout_matrix = matrix(c(1,4,
                         1,4,
                         2,4,
                         2,4,
                         3,4,
                         3,5),nrow=6,byrow=T)
layout(layout_matrix, heights=c(1,1,1,1,.6,1.4))


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

write(paste('Combined P(G) for Launched T-I: ',percent(combined_ti_forest$mean[5],digits=1),'\n',sep=''),text_stats_path,append=T)

write(paste('Combined P(G) for all active clinical T-I: ',
            percent(sum(active_ti_forest$numerator[active_ti_forest$label %in% active_clinical$cat])/sum(active_ti_forest$denominator[active_ti_forest$label %in% active_clinical$cat]),digits=1),
            ' (',sum(active_ti_forest$numerator[active_ti_forest$label %in% active_clinical$cat]),
            '/',sum(active_ti_forest$denominator[active_ti_forest$label %in% active_clinical$cat]),')','\n',sep=''),text_stats_path,append=T)

write(paste('Combined P(G) for all historical clinical T-I: ',
            percent(sum(hist_ti_forest$numerator[hist_ti_forest$label %in% active_clinical$cat])/sum(hist_ti_forest$denominator[hist_ti_forest$label %in% active_clinical$cat]),digits=1),
            ' (',sum(hist_ti_forest$numerator[hist_ti_forest$label %in% active_clinical$cat]),
            '/',sum(hist_ti_forest$denominator[hist_ti_forest$label %in% active_clinical$cat]),')','\n',sep=''),text_stats_path,append=T)

rbind(cbind(hist_ti_forest, mode='historical'),
      cbind(active_ti_forest, mode='active'),
      cbind(combined_ti_forest, mode='combined')) %>%
  relocate(mode) %>% 
  filter(label != 'Historical') %>%
  ungroup() %>%
  select(-num) %>%
  rename(phase=label, supported=numerator, total=denominator, pg=mean, pg_l95=l95, pg_u95=u95) %>%
  select(mode, phase, supported, total, pg, pg_l95, pg_u95) -> pg_phase_mode
write_supp_table(pg_phase_mode, "Proportion of target-indication pairs with genetic support by phase and active/historical status.")

#### 1B - by association source
plot_forest(assoc_source_rr_forest, xlims=c(0,5), xstyle='ratio', mar=c(2.5, 8, 3, 8))
mtext(side=1, line=1.6, text='RS')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

assoc_source_rr_forest %>%
  select(source=label, rs=mean, rs_l95 = l95, rs_u95 = u95, approved=numerator, supported=denominator) -> assoc_source_rr_out
write_supp_table(assoc_source_rr_out, "Relative success by source of germline genetic evidence.")

##### 1C - by threshold
roc_col = '#FFAA00'
omim_col = '#00CDCD'
xlims = c(1000, 0)
ylims = c(0.8, 4.0)
par(mar=c(4,4,3,6))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=2, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=0:50/10, labels=NA, tck=-0.025)
axis(side=2, at=0:5, tck=-0.05, las=2)
abline(h=1, lty=2)
axis(side=1, at=0:10*100, labels=formatC(0:10*100,big.mark=','))
roc_l2g %>% filter(!is.na(roc_l2g$l95)) -> roc_plot
points(y=roc_plot$mean, x=roc_plot$gensup_all, type='l', lwd=2, col=roc_col)
points(y=roc_plot$mean, x=roc_plot$gensup_all, pch=20, cex=0.75, col=roc_col)
polygon(y=c(roc_plot$l95,rev(roc_plot$u95)), x=c(roc_plot$gensup_all, rev(roc_plot$gensup_all)), col=alpha(roc_col,ci_alpha), border=roc_col, lwd=0.5, lty=3)
to_label = roc_plot$share %in% c(0.25, 0.5, 0.75)
points(y=roc_plot$mean[to_label], x=roc_plot$gensup_all[to_label], pch=20, cex=0.25, col='black')
text(y=roc_plot$mean[to_label], x=roc_plot$gensup_all[to_label], pos=3, labels=roc_plot$share[to_label])
mtext(side=2, line=2, text='RS')
mtext(side=1, line=2.5, text='N genetically supported T-I pairs')
mtext(side=4, at=roc_plot$mean[1], text='OTG', las=2, line=0.25, col=roc_col)
abline(h=assoc_source_rr_forest$mean[assoc_source_rr_forest$label=='OMIM'], col=omim_col, lwd=1)
abline(h=assoc_source_rr_forest$l95[assoc_source_rr_forest$label=='OMIM'], col=omim_col, lwd=0.5, lty=3)
mtext(side=4, at=assoc_source_rr_forest$mean[assoc_source_rr_forest$label=='OMIM'], text='OMIM', col=omim_col, las=2, line=0.25, cex = 0.75)
mtext(side=4, at=assoc_source_rr_forest$l95[assoc_source_rr_forest$label=='OMIM']-.2, text='95%CI', col=omim_col, las=2, line=0.25, cex = 0.75)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

roc_l2g %>%
  select(-numerator, -denominator) -> roc_l2g_out
write_supp_table(roc_l2g_out, "Relative success for Open Targets Genetics associations as a function of locus to gene (L2G) share.")

master_forest$y = nrow(master_forest):1
master_forest$subpanel = c(rep('Year',sum(year_rrs$abbr=='c')),
                           rep('Gene count',nrow(gc_rrs)),
                           rep('Beta',nrow(beta_rrs)),
                           rep('Odds ratio',nrow(or_rrs)),
                           rep('MAF',nrow(maf_rrs)))
master_forest %>% relocate(subpanel) -> master_forest

write_supp_table(master_forest, "Relative success for GWAS Catalog associations as a function of year of discovery, gene count, beta, odds ratio, and minor allele frequency.")

yearfirst_logit_data %>%
  select(ti_uid, gene, indication_mesh_id, indication_mesh_term, ccat, similarity, assoc_year, assoc_source, assoc_mesh_id, assoc_mesh_term, original_trait, original_link, l2g_share) -> yearfirst_logit_data_out
write_supp_table(yearfirst_logit_data, 'Launched status versus year of discovery for GWAS-supported target-indication pairs (data for logit model).')

master_forest$label[grepl('^All',master_forest$label)] = 'All'
master_forest$label = gsub('^(beta |OR |MAF )','',master_forest$label)
master_forest %>%
  group_by(subpanel) %>%
  summarize(.groups='keep',
            miny = min(y),
            maxy = max(y),
            midy = mean(y)) -> master_tranches

tranche_offset = 7

plot_forest(master_forest, xlims=c(0,4), xstyle='ratio', mar=c(3,8,3,6))
for (i in 1:nrow(master_tranches)) {
  axis(side=2, line=tranche_offset, at=c(master_tranches$miny[i]-.25, master_tranches$maxy[i]+.25), tck=0.05, labels=NA)
  mtext(side=2, at=master_tranches$midy[i], text=master_tranches$subpanel[i], line=tranche_offset+.25)
}
abline(h=master_tranches$maxy+0.5,lwd=0.5)
mtext(side=1, line=1.6, text='RS')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1




pp %>%
  group_by(indication_mesh_id) %>%
  summarize(.groups='keep',
            n_human_target = length(unique(na.omit(gene)))) %>% 
  ungroup() -> targets_per_indic

read_tsv('data/drugs_per_indic.tsv', col_types=cols()) -> drugs_per_indic 

indic %>%
  filter(!is.na(indication_mesh_term)) %>%
  select(indication_mesh_id, indication_mesh_term) %>%
  left_join(sim, by=c('indication_mesh_id' = 'meshcode_a')) %>%
  filter(comb_norm >= 0.8) %>%
  left_join(assoc, by=c('meshcode_b' = 'mesh_id')) %>%
  filter(l2g_share >= 0.5 | source != 'OTG') %>%
  rename(assoc_mesh_id=meshcode_b, assoc_mesh_term=mesh_term) %>%
  group_by(indication_mesh_id, indication_mesh_term, assoc_mesh_id, assoc_mesh_term, gene) %>%
  summarize(.groups='keep', n_entries=n()) %>%
  ungroup() %>%
  group_by(indication_mesh_id) %>%
  summarize(.groups='keep', 
            n_supported_genes=length(unique(na.omit(gene)))) %>%
  ungroup() -> assocs_per_indic

targets_per_indic %>%
  left_join(drugs_per_indic, by=c('indication_mesh_id'='mesh_id')) %>%
  left_join(assocs_per_indic, by=c('indication_mesh_id')) %>%
  mutate(n_supported_genes = replace_na(n_supported_genes, 0)) %>%
  mutate(logbin = ceil(log10(n_supported_genes))) %>%
  mutate(logbin = replace(logbin, logbin==-Inf, -1)) %>%
  filter(!is.na(indication_mesh_id)) -> indic_stats
indic_stats$indication_mesh_term = mesh_best_names$labeltext[match(indic_stats$indication_mesh_id, mesh_best_names$id)]

indic_stats %>%
  filter(n_human_target > 0) %>%
  group_by(logbin) %>%
  summarize(.groups='keep', n_indic=n()) %>%
  ungroup() %>%
  mutate(bin_name = case_when(logbin == -1 ~ '0',
                              logbin == 0 ~ '1',
                              logbin == 1 ~ '2-10',
                              logbin == 2 ~ '11-100',
                              logbin == 3 ~ '101-999',
                              logbin == 4 ~ '1,000+')) %>%
  arrange(logbin) -> gene_bins

indic_stats %>%
  filter(n_human_target > 0) %>%
  filter(n_supported_genes == 0) %>%
  arrange(desc(n_drugs)) -> indications_without_supported_genes

par(mar=c(3,3,3,1))
xlims = c(-1.5, 4.5)
ylims = c(0, 450)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=0:5*100, labels=NA, tck=-0.02)
axis(side=2, at=0:5*100, las=2, lwd=0, lwd.ticks=0, line=-0.5, cex.axis=0.9)
axis(side=1, at=xlims, lwd=1, lwd.ticks=0, labels=NA)
mtext(side=1, line=0.25, at=gene_bins$logbin, text=gene_bins$bin_name, cex=0.7)
#axis(side=1, at=gene_bins$logbin, labels=gene_bins$bin_name, lwd.ticks=0, lwd=0, line=-0.5)
barwidth=0.4
rect(xleft=gene_bins$logbin-barwidth, xright=gene_bins$logbin+barwidth, ybottom=rep(0,nrow(gene_bins)), ytop=gene_bins$n_indic, border=NA, col='#A7A7A7')
mtext(side=1, line=1.6, text='N genes associated with similar traits', cex=1)
mtext(side=2, line=2.0, text='N indications', cex=1)
mtext(letters[panel], side=3, cex=2, adj = 0.00, line = 0.5)
panel = panel + 1

write(paste('Indications with 0 supported genes: ',gene_bins$n_indic[gene_bins$logbin==-1],'\n',sep=''),text_stats_path,append=T)

write_supp_table(gene_bins, "Number of indications by number of supported genes.")

write_supp_table(indications_without_supported_genes, "Indications developed in Pharmaprojects for which there exist no genetically supported targets.")

unecessary_message = dev.off()














####
# Additional Figure S2/3 staging
####

cat(file=stderr(), 'done.\nCalculating oncology, orphan, and ROC sim for Figure S2...')

intogen_genes = read.table('data/intogen_genes.tsv',sep='\t',header=T)

onco_rrs = tibble(y=16:1,
                  label=c('All indications','oncology','    oncogenes','    tumor suppressors','    unknown','non-oncology','oncology','    oncogenes','    tumor suppressors','    unknown', 'non-oncology','oncology','    oncogenes','    tumor suppressors','    unknown','non-oncology'),
                  areas = c('All',rep('C04',15)),
                  area_filter = c('All',rep(c(rep('only',4),'non'),3)),
                  assoc_source = c('intOGen',rep(c('intOGen','OMIM','GWAS'),each=5)),
                  intogen_mechanism = c('All',rep(c('All','oncogene','tumor suppressor','unknown','All'),3)))
onco_rrs %>%
  group_by(assoc_source) %>%
  summarize(.groups='keep', miny=min(y), maxy=max(y), midy=mean(y)) %>%
  arrange(desc(maxy)) -> y_tranches

for (i in 1:nrow(onco_rrs)) {
  if (onco_rrs$assoc_source[i] == 'intOGen') {
    orig_obj = hist_ti_intogen
  } else if ( onco_rrs$assoc_source[i] == 'OMIM') {
    orig_obj = hist_ti_omim
  } else if (onco_rrs$assoc_source[i] == 'GWAS') {
    orig_obj = hist_ti_gwas
  }
  area_subset = subset_by_area(orig_obj, topl=onco_rrs$areas[i], filter=onco_rrs$area_filter[i])
  if (onco_rrs$intogen_mechanism[i]=='All') {
    gene_subset = area_subset
  } else {
    gene_subset = area_subset[area_subset$gene %in% intogen_genes$gene[intogen_genes$mechanism==onco_rrs$intogen_mechanism[i]],]
  }
  rr_obj = advancement_rr(gene_subset)
  onco_rrs[i,c('mean','l95','u95')] = rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  onco_rrs[i,c('numerator')]        = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  onco_rrs[i,c('denominator')]      = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
}

roc_col = '#FFAA00'
roc_sim = tibble(thresh=seq(0.1, 1.0, 0.05))
for (i in 1:nrow(roc_sim)) {
  rr_obj = advancement_rr(hist_ti, threshold=roc_sim$thresh[i])
  roc_sim[i,c('mean','l95','u95')] = rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  roc_sim[i,c('numerator')]        = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  roc_sim[i,c('denominator')]      = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
  roc_sim[i,c('gensup_all')]       = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
  roc_sim[i,c('gensup_launched')]  = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
}

# hist_ti_firstyear = pipeline_best(merge2, basis='ti', phase='historical', associations='OTG', otg_subcat='GWAS Catalog', firstyear=T, minusomim=T)
# hist_ti_firstyear %>%
#   filter(ccat=='Launched') -> hist_ti_firstyear_launched




#### 

hist_ti_omim_orphan = subset_by_area(hist_ti_omim, topl='all', orphan='yes')
hist_ti_omim_nonorphan = subset_by_area(hist_ti_omim, topl='all', orphan='no')
hist_ti_otg_orphan = subset_by_area(hist_ti_otg, topl='all', orphan='yes')
hist_ti_otg_nonorphan = subset_by_area(hist_ti_otg, topl='all', orphan='no')
orphan_forest = tibble(source = rep(c('OMIM','OTG'),each=2),
                       orphan = rep(c('orphan','non-orphan'),2),
                                pipeline_obj=c('hist_ti_omim_orphan','hist_ti_omim_nonorphan','hist_ti_otg_orphan','hist_ti_otg_nonorphan')) %>%
  mutate(label = paste0(source,' ',orphan)) %>%
  mutate(y=max(row_number()) - row_number() + 1)
for (i in 1:nrow(orphan_forest)) {
  pipeline_obj = get(orphan_forest$pipeline_obj[i])
  rr_obj = advancement_rr(pipeline_obj)
  orphan_forest[i,c('mean','l95','u95')] = rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  orphan_forest[i,c('numerator')]        = rr_obj[rr_obj$phase=='I-Launch',c('x_yes')]
  orphan_forest[i,c('denominator')]      = rr_obj[rr_obj$phase=='I-Launch',c('n_yes')]
}

cat(file=stderr(), 'done.\nCreating Figure S2...')

resx=300
png(paste0(output_path,'/figure-s2.png'),width=6.5*resx,height=8*resx,res=resx)

layout_matrix = matrix(c(1,2,
                         3,4,
                         3,5,
                         6,7),nrow=4,byrow=T)
layout(layout_matrix, heights=c(.8,.6, .6, 1.2))

panel = 1

xlims = c(max(roc_sim$gensup_all)+500, 0)
ylims = c(0, 3.5)
par(mar=c(3,4,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=2, at=0:40/10, labels=NA, tck=-0.025,las=2)
axis(side=2, at=0:4, tck=-0.05,las=2)
abline(h=1, lty=2)
axis(side=1, at=0:15*1000, labels=NA, tck=-0.025) 
axis(side=1, at=0:3*5000, labels=NA, tck=-0.05)
axis(side=1, at=0:3*5000, line=-0.5, lwd=0, labels=formatC(0:3*5000,big.mark=',',format='d'))
points(y=roc_sim$mean, x=roc_sim$gensup_all, type='l', lwd=2, col=roc_col)
points(y=roc_sim$mean, x=roc_sim$gensup_all, pch=20, cex=0.75, col=roc_col)
polygon(y=c(roc_sim$l95,rev(roc_sim$u95)), x=c(roc_sim$gensup_all, rev(roc_sim$gensup_all)), col=alpha(roc_col,ci_alpha), border=NA)
to_label = roc_sim$thresh %in% c(0.2, 0.4, 0.6, 0.8)
text(y=roc_sim$mean[to_label], x=roc_sim$gensup_all[to_label], pos=3, labels=roc_sim$thresh[to_label])
points(y=roc_sim$mean[to_label], x=roc_sim$gensup_all[to_label], pch=1)
mtext(side=2, line=1.5, text='RS', cex=0.75)
mtext(side=1, line=1.6, text='N genetically supported T-I pairs', cex=0.75)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(roc_sim, "Relative success as a function of indication-association similarity threshold.")

plot_forest(orphan_forest, xlims=c(0,4), xstyle='ratio', mar=c(4,12,4,6))
mtext(side=1, line=2.0, text='RS')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(orphan_forest, "Relative success for OMIM and OTG associations by orphan status.")

plot_forest(onco_rrs, xlims=c(0,5), xstyle='ratio', mar=c(3,12,3,6))
mtext(side=1, line=2.0, text='RS')
overhang = 0.35
tranche_line = 8.75
for (i in 1:nrow(y_tranches)) {
  axis(side=2, line=tranche_line, at=c(y_tranches$maxy[i]+overhang,y_tranches$miny[i]-overhang), labels=NA, tck=0.025)
  axis(side=2, line=tranche_line, at=y_tranches$midy[i], labels=y_tranches$assoc_source[i], lwd=0, cex.axis=1.2)
}
abline(h=y_tranches$maxy+0.5, lwd=0.5)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(onco_rrs, "Relative success for somatic vs. germline support in oncology.")

plot_forest(year_rrs[year_rrs$abbr=='a',], xlims=c(0,4), xstyle='ratio', mar=c(3,6,3,6))
mtext(side=1, line=1.6, text='RS')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(year_rrs[year_rrs$abbr=='a',], "Relative success for GWAS Catalog associations by year of discovery, without removing replications or OMIM.")

plot_forest(year_rrs[year_rrs$abbr=='b',], xlims=c(0,4), xstyle='ratio', mar=c(3,6,3,6))
mtext(side=1, line=1.6, text='RS')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(year_rrs[year_rrs$abbr=='b',], "Relative success for GWAS Catalog associations by year, removing replications but not removing OMIM.")

yearfirst_scatter_pbest = pipeline_best(merge2, 
                                    phase='historical', 
                                    basis='ti', 
                                    associations='OTG', # note no otg_subcat restriction here
                                    min_year=2007,
                                    max_year=2022,
                                    firstyear=T,
                                    minusomim=T, verbose=F)
yearfirst_scatter_data = subset(yearfirst_scatter_pbest, target_status=='genetically supported target' & cat %in% c('Launched',active_clinical$cat))


set.seed(1)
yearfirst_scatter_data %>%
  filter(cat=='Launched' & !is.na(year_launch)) %>%
  select(ti_uid, 
         gene, 
         min_assoc_year,
         year_launch) %>%
  group_by(min_assoc_year, year_launch) %>%
  mutate(n_points = n()) %>%
  ungroup() %>%
  mutate(xjit = case_when(n_points > 1 ~ jitter(min_assoc_year,amount=.125),
                          n_points==1 ~ min_assoc_year)) -> yearfirst_scatter

yearfirst_scatter %>%
  arrange(gene, year_launch) %>%
  group_by(gene) %>%
  slice(1) %>%
  ungroup() %>%
  filter(year_launch - 5 >= min_assoc_year) %>%
  select(gene, min_assoc_year, year_launch) %>%
  mutate(pos=2) %>%
  mutate(pos = case_when(gene=='ERBB4' ~ 4,
                         gene=='BCL2' ~ 4,
                         gene=='CASR' ~ 1,
                         TRUE ~ 2)) -> yearfirst_highlights

par(mar=c(4,4.5,3,1))
xlims = c(2000, 2023)
ylims = c(2000, 2023)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
abline(a=0,b=1,col='gray')
axis(side=1, at=xlims, lwd=1, lwd.ticks=0, labels=NA)
axis(side=1, lwd=0, lwd.ticks=1, at=seq(2000,2023,1), labels=NA, tck=-0.02)
axis(side=1, lwd=0, lwd.ticks=1, at=seq(2000,2023,5), labels=NA, tck=-0.05)
axis(side=1, lwd=0, lwd.ticks=0, at=seq(2000,2023,5), line=-0.25)
axis(side=2, at=xlims, lwd=1, lwd.ticks=0, labels=NA)
axis(side=2, lwd=0, lwd.ticks=1, at=seq(2000,2023,1), labels=NA, tck=-0.02)
axis(side=2, lwd=0, lwd.ticks=1, at=seq(2000,2023,5), labels=NA, tck=-0.05)
axis(side=2, lwd=0, lwd.ticks=0, at=seq(2000,2023,5), line=-0.25, las=2)
points(x=yearfirst_scatter$xjit, y=yearfirst_scatter$year_launch, pch=20)
mtext(side=1, line=2.5, text='Association year')
mtext(side=2, line=3, text='Launch year')
par(xpd=T)
points(x=yearfirst_highlights$min_assoc_year, yearfirst_highlights$year_launch, pch=1, cex=1.5)
text(x=yearfirst_highlights$min_assoc_year, yearfirst_highlights$year_launch, labels=yearfirst_highlights$gene, pos=yearfirst_highlights$pos, font=3)
par(xpd=F)

yearfirst_scatter %>%
  select(-n_points, -xjit) -> yearfirst_scatter_output
write_supp_table(yearfirst_scatter_output, "Years of association and launch for genetically supported launched target-indication pairs.")

mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


par(mar=c(6,4,6,1))
xlims = c(2007, 2022)
ylims = c(0, 1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=2005:2022, las=2)
mtext(side=1, line=3, text='First year of genetic association')
axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2)
mtext(side=2, line=3, text='P(Launched|Supported)')
points(yearfirst_bins$assoc_year, yearfirst_bins$p_launched, type='l', lwd=3)
polygon(x=c(yearfirst_bins$assoc_year, rev(yearfirst_bins$assoc_year)), y=c(yearfirst_bins$l95, rev(yearfirst_bins$u95)), col='#A9A9A988', border=NA)
mtext(side=3, at=yearfirst_bins$assoc_year, text=paste0(yearfirst_bins$n_launched,'/',yearfirst_bins$n_i_l), las=2, line=0.5)
mtext(letters[panel], side=3, cex=2, adj = -0.15, line = 1.5)
panel = panel + 1

write_supp_table(yearfirst_bins, "Proportion of supported target-indication pairs that are launched, by year of discovery.")

unecessary_message = dev.off()







resx=300
png(paste0(output_path,'/figure-s3.png'),width=6.5*resx,height=4.5*resx,res=resx)

# all synonyms for all MeSH terms
all_vocab = read_tsv('data/mesh_all_vocab.tsv.gz', col_types=cols())
# prepare to match them to Nelson 2015
all_vocab %>%
  mutate(labeltext = tolower(labeltext)) %>%
  group_by(labeltext) %>%
  slice(1) %>%
  ungroup() -> vocab_match


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

unnecessary_message = dev.off()

# 
# # examples of what is in signs & symptoms
# pp %>% inner_join(indic_topl_match, by=c('indication_mesh_id','indication_mesh_term')) %>% filter(topl=='C23') %>%
#   group_by(indication_mesh_id, indication_mesh_term) %>%
#   summarize(.groups='keep', n=n()) %>%
#   arrange(desc(n))




cat(file=stderr(), 'done.\nCreating Figure S4...')

resx=300
png(paste0(output_path,'/figure-s4.png'),width=6.5*resx,height=3.5*resx,res=resx)

# read data
t2d_omim = read_tsv('data/t2d/omim_t2d.tsv', col_types=cols())
t2d_gwas = read_tsv('data/t2d/t2d_gwas_vujkovic_suzuki.tsv', col_types=cols())
pp_dia = read_tsv('data/t2d/pp_diabetes.tsv', col_types=cols())
curated = read_tsv('data/t2d/drug_gene_curation_results.tsv', col_types=cols())

t2d_omim_out = t2d_omim %>% arrange(gene, mesh_id)
t2d_gwas_out = t2d_gwas %>% arrange(gene)
write_supp_table(t2d_omim_out, "T2D genes with OMIM support.")
write_supp_table(t2d_gwas_out, "T2D genes with GWAS support by novelty status.")

pp_dia$keep = curated$keep[match(pp_dia$gene, curated$gene)]
pp_dia$keep[is.na(pp_dia$keep)] = TRUE
pp_dia$drug = curated$drug[match(pp_dia$gene, curated$gene)]
pp_dia$comments = curated$comments[match(pp_dia$gene, curated$gene)]
pp_dia$gensup_omim = pp_dia$gene %in% t2d_omim$gene
pp_dia$gensup_estab = pp_dia$gene %in% t2d_gwas$gene[!t2d_gwas$novel]
pp_dia$gensup_novel = pp_dia$gene %in% t2d_gwas$gene[t2d_gwas$novel]
pp_dia$support = 'none'
pp_dia$support[pp_dia$gensup_omim] = 'omim'
pp_dia$support[pp_dia$gensup_novel & !pp_dia$gensup_omim] = 'novel gwas'
pp_dia$support[pp_dia$gensup_estab & !pp_dia$gensup_omim] = 'established gwas' 
pp_dia$support[pp_dia$gensup_novel & pp_dia$gensup_omim] = 'omim, novel gwas'
pp_dia$support[pp_dia$gensup_estab & pp_dia$gensup_omim] = 'omim, established gwas' 
pp_dia$gwas_source = t2d_gwas$source[match(pp_dia$gene, t2d_gwas$gene)]
pp_dia$gwas_source[is.na(pp_dia$gwas_source)] = 'None'


pp_dia %>%
  filter(keep) %>%
  select(gene, indication_mesh_id, indication_mesh_term, maxphasenum, maxphase, drug, support) -> pp_dia_output

write_supp_table(pp_dia_output, "T2D drug development programs by genetic support status.")

pp_dia %>%
  filter(keep) %>%
  group_by(maxphasenum, maxphase) %>%
  summarize(.groups='keep',
            n_omim = sum(gensup_omim),
            n_estab = sum(gensup_estab),
            n_novel = sum(gensup_novel),
            n_all = sum(gensup_omim | gensup_estab | gensup_novel),
            n_total = n()) %>%
  ungroup() %>%
  arrange(maxphasenum, maxphase) -> forest

forest$y = max(forest$maxphasenum) - forest$maxphasenum + 1
forest[,c('omim_mean','omim_l95','omim_u95')] = binom.confint(x=forest$n_omim, n=forest$n_total, method='wilson')[,c('mean','lower','upper')]
forest[,c('novel_mean','novel_l95','novel_u95')] = binom.confint(x=forest$n_novel, n=forest$n_total, method='wilson')[,c('mean','lower','upper')]
forest[,c('estab_mean','estab_l95','estab_u95')] = binom.confint(x=forest$n_estab, n=forest$n_total, method='wilson')[,c('mean','lower','upper')]
forest[,c('all_mean','all_l95','all_u95')] = binom.confint(x=forest$n_all, n=forest$n_total, method='wilson')[,c('mean','lower','upper')]

write_supp_table(forest, "T2D proportion of programs with genetic support by development phase.")

panel = 1

par(mar=c(2,0.5,4,4))
layout_matrix = matrix(1:6,nrow=2,byrow=T)
layout(layout_matrix,widths=c(.4,1,1))
ylims = range(forest$y) + c(-0.5, 0.5)
xlims = c(0,.6)
panel = 1
# legend panel
plot(NA,NA,xlim=xlims,ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
mtext(side=4,adj=1,line=3,at=forest$y,text=forest$maxphase,las=2)
# omim
plot(NA,NA,xlim=xlims,ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=1, at=0:100/100, labels=NA, tck=-0.025)
axis(side=1, at=0:10/10, labels=NA, tck=-0.05)
axis(side=1, at=0:3/5, labels=percent(0:3/5),lwd=0,line=-0.5)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
points(x=forest$omim_mean, y=forest$y, pch=19)
segments(x0=forest$omim_l95, x1=forest$omim_u95, y0=forest$y, lwd=1.5)
mtext(side=4, at=forest$y, line=0.25, text=paste(forest$n_omim,forest$n_total,sep='/'),las=2,cex=0.85)
mtext(side=3, line=0, text='OMIM')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.25)
panel = panel + 1
# gwas estab
plot(NA,NA,xlim=xlims,ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=1, at=0:100/100, labels=NA, tck=-0.025)
axis(side=1, at=0:10/10, labels=NA, tck=-0.05)
axis(side=1, at=0:3/5, labels=percent(0:3/5),lwd=0,line=-0.5)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
points(x=forest$estab_mean, y=forest$y, pch=19)
segments(x0=forest$estab_l95, x1=forest$estab_u95, y0=forest$y, lwd=1.5)
mtext(side=4, at=forest$y, line=0.25, text=paste(forest$n_estab,forest$n_total,sep='/'),las=2,cex=0.85)
mtext(side=3, line=0, text='GWAS established')
par(xpd=T)
mtext(side=1, line=1.75, at=max(ylims), text='Proportion of drug mechanisms with human genetic support')
par(xpd=F)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.25)
panel = panel + 1
par(mar=c(4,0.5,2,4))
# legend panel
plot(NA,NA,xlim=xlims,ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
mtext(side=4,adj=1,line=3,at=forest$y,text=forest$maxphase,las=2)
# gwas novel
plot(NA,NA,xlim=xlims,ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=1, at=0:100/100, labels=NA, tck=-0.025)
axis(side=1, at=0:10/10, labels=NA, tck=-0.05)
axis(side=1, at=0:3/5, labels=percent(0:3/5),lwd=0,line=-0.5)
mtext(side=1, line=2, at=max(xlims)+.1, text='Proportion of drug mechanisms with human genetic support', cex=1)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
points(x=forest$novel_mean, y=forest$y, pch=19)
segments(x0=forest$novel_l95, x1=forest$novel_u95, y0=forest$y, lwd=1.5)
mtext(side=4, at=forest$y, line=0.25, text=paste(forest$n_novel,forest$n_total,sep='/'),las=2,cex=0.85)
mtext(side=3, line=0, text='GWAS novel')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.25)
panel = panel + 1
plot(NA,NA,xlim=xlims,ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=1, at=0:100/100, labels=NA, tck=-0.025)
axis(side=1, at=0:10/10, labels=NA, tck=-0.05)
axis(side=1, at=0:3/5, labels=percent(0:3/5),lwd=0,line=-0.5)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
points(x=forest$all_mean, y=forest$y, pch=19)
segments(x0=forest$all_l95, x1=forest$all_u95, y0=forest$y, lwd=1.5)
mtext(side=4, at=forest$y, line=0.25, text=paste(forest$n_all,forest$n_total,sep='/'),las=2,cex=0.85)
mtext(side=3, line=0, text='All sources')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.25)
panel = panel + 1


unnecessary_message = dev.off()













cat(file=stderr(), 'done.\nCreating Figure 2...')

read_tsv('data/areas.tsv', col_types=cols()) %>%
  mutate(filter='only') %>%
  select(topl, area, color, filter) %>%
  add_row(topl='ALL',area='all',color='#000000',filter='', .before=1) -> areas_all

if(exists('rr_table_areas')) {
  rm(rr_table_areas)
}
for (i in 1:nrow(areas_all)) {
  hist_ti_area = subset_by_area(hist_ti, areas_all$topl[i], areas_all$filter[i])
  area_rr_obj = advancement_rr(hist_ti_area)
  areas_all[i,c('rs_mean','rs_l95','rs_u95')]         = area_rr_obj[area_rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  areas_all[i,'n_hist_ti'] = sum(area_rr_obj[area_rr_obj$phase=='I-Launch',c('n_yes','n_no')])
  
  # discuss whether we still need all of these extra columns and separate rr_table_areas table
  areas_all[i,c('nosup_clinical','nosup_launched')]   = area_rr_obj[area_rr_obj$phase=='I-Launch',c('n_no','x_no')]
  areas_all[i,c('gensup_clinical','gensup_launched')] = area_rr_obj[area_rr_obj$phase=='I-Launch',c('n_yes','x_yes')]
  areas_all[i,c('rs_i_l_mean','rs_i_l_l95','rs_i_l_u95')]          = area_rr_obj[area_rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')]
  areas_all[i,c('rs_i_ii_mean','rs_i_ii_l95','rs_i_ii_u95')]       = area_rr_obj[area_rr_obj$phase=='I',c('rs_mean','rs_l','rs_u')]
  areas_all[i,c('rs_ii_iii_mean','rs_ii_iii_l95','rs_ii_iii_u95')] = area_rr_obj[area_rr_obj$phase=='II',c('rs_mean','rs_l','rs_u')]
  areas_all[i,c('rs_iii_a_mean','rs_iii_a_l95','rs_iii_a_u95')]    = area_rr_obj[area_rr_obj$phase=='III',c('rs_mean','rs_l','rs_u')]
  
  area_rr_obj$area   = areas_all$area[i]
  area_rr_obj$topl   = areas_all$topl[i]
  area_rr_obj$filter = areas_all$filter[i]
  if(exists('rr_table_areas')) {
    rr_table_areas = rbind(rr_table_areas, area_rr_obj)
  } else {
    rr_table_areas = area_rr_obj
  }
}

areas_all %>%
  filter(!is.na(rs_mean)) %>%
  arrange(filter, desc(rs_mean)) %>%
  mutate(y = max(row_number()) - row_number() + 1) -> areas_all

rr_table_areas$y = areas_all$y[match(rr_table_areas$area, areas_all$area)]
rr_table_areas$color = areas_all$color[match(rr_table_areas$area, areas_all$area)]

areas_all %>%
  filter(area != 'all') -> areas

rr_table_areas %>%
  select(area, topl, phase:fraction) -> rr_table_output

write_supp_table(rr_table_output, "Relative success by therapy area and phase transition.")

# and prep some stats for text and supplement
cmh_array = array(data=rep(0,2*2*(nrow(areas))), dim = c(2,2,(nrow(areas))))
for (i in 1:nrow(areas)) {
  cmh_array[1,1,i] = as.numeric(areas[i,'nosup_clinical'])
  cmh_array[1,2,i] = as.numeric(areas[i,'nosup_launched'])
  cmh_array[2,1,i] = as.numeric(areas[i,'gensup_clinical'])
  cmh_array[2,2,i] = as.numeric(areas[i,'gensup_launched'])
}
cmh_result = suppressMessages(cmh.test(cmh_array))
extract_cmh_p = function(cmh_obj) {
  cmh_minp = 1e-15 # 1 - pchisq(q=70, df=1) = 1.1e-16, this is the lowest non-zero p value it ever returns, after that all 0. so we can just say P < 1e-15
  chisq_mh = cmh_obj$parameter['CMH statistic']
  df = cmh_obj$parameter['df']
  p = 1 - pchisq(q=chisq_mh, df=df)
  return (max(p,cmh_minp))
}

write(paste('CMH test on I-Launched RS values across ',dim(cmh_array)[3],' therapy areas: P = ',formatC(extract_cmh_p(cmh_result),format='e',digits=1),'\n',sep=''),text_stats_path,append=T)

p_s_data = data.frame(n_success=areas$nosup_launched+areas$gensup_launched, n_failure=areas$nosup_clinical + areas$gensup_clinical)
p_s_mat = as.matrix(p_s_data, colnames=c('success','failure'))
p_s_chisq_p = suppressWarnings(chisq.test(p_s_mat)$p.value)
p_s_bconf_obj = binom.confint(x=areas$nosup_launched+areas$gensup_launched, n=areas$nosup_launched+areas$gensup_launched+areas$nosup_clinical+areas$gensup_clinical, method='wilson')
areas[,c('p_s_mean','p_s_l95','p_s_u95')] = p_s_bconf_obj[,c('mean','lower','upper')]
p_g_data = data.frame(n_gensup=areas$gensup_launched + areas$gensup_clinical, n_nosup=areas$nosup_launched + areas$nosup_clinical)
p_g_mat = as.matrix(p_g_data, colnames=c('gensup','nosup'))
p_g_chisq_p = suppressWarnings(chisq.test(p_g_mat)$p.value)
p_g_bconf_obj = binom.confint(x=areas$gensup_clinical+areas$gensup_launched, n=areas$nosup_launched+areas$gensup_launched+areas$nosup_clinical+areas$gensup_clinical, method='wilson')
areas[,c('p_g_mean','p_g_l95','p_g_u95')] = p_g_bconf_obj[,c('mean','lower','upper')]

indic %>% filter(genetic_insight != 'none') -> indic_insight





# here we go back to pp table in order to get all indications regardless of insight
pp %>%
  filter(ccat != 'Preclinical' & gene != '' & indication_mesh_id != '') %>%
  group_by(gene) %>%
  summarize(.groups='keep', 
            n_clinical_indic = length(unique(indication_mesh_id)),
            n_launched_indic = length(unique(indication_mesh_id[ccat=='Launched']))) -> ti_stats
pp %>%
  filter(hcat=='Launched' & gene != '' & indication_mesh_id != '') %>%
  distinct(gene, indication_mesh_id) -> ti_launched
ti_launched %>%
  select(gene, meshcode_a = indication_mesh_id) -> a1
ti_launched %>%
  select(gene, meshcode_b = indication_mesh_id) -> a2
a1 %>%
  inner_join(a2, by=c('gene'='gene')) %>%
  filter(meshcode_a < meshcode_b) -> ti_simpairs
ti_simpairs %>%
  left_join(sim, by=c('meshcode_a' = 'meshcode_a', 'meshcode_b' = 'meshcode_b')) -> ti_sim
ti_sim$comb_norm[is.na(ti_sim$comb_norm)] = 0 # in current dataset this never arises - occurs only with corrupted mesh codes, which we fixed
ti_launched %>%
  group_by(gene) %>%
  summarize(.groups='keep', 
            n_launched_indic = length(unique(indication_mesh_id))) -> target_stats
ti_sim %>%
  group_by(gene) %>%
  summarize(.groups='keep', 
            meansim=mean(comb_norm)) -> sim_stats

target_stats$meansim = sim_stats$meansim[match(target_stats$gene, sim_stats$gene)]
target_stats$meansim[target_stats$n_launched_indic==1] = 1.00




# main version with all programs all time
pp_launched_alltime = read_tsv('data/pp_launched_alltime.tsv', col_types=cols())
pp_launched_alltime %>%
  select(gene, meshcode_a = indication_mesh_id) -> a1
pp_launched_alltime %>%
  select(gene, meshcode_b = indication_mesh_id) -> a2
a1 %>%
  inner_join(a2, by=c('gene'='gene')) %>%
  filter(meshcode_a < meshcode_b) -> ti_simpairs
ti_simpairs %>%
  left_join(sim, by=c('meshcode_a' = 'meshcode_a', 'meshcode_b' = 'meshcode_b')) -> ti_sim
ti_sim$comb_norm[is.na(ti_sim$comb_norm)] = 0 # in current dataset this never arises - occurs only with corrupted mesh codes, which we fixed
pp_launched_alltime %>%
  group_by(gene) %>%
  summarize(.groups='keep', 
            n_launched_indic = length(unique(indication_mesh_id))) %>%
  ungroup() -> target_stats
ti_sim %>%
  group_by(gene) %>%
  summarize(.groups='keep', 
            meansim=mean(comb_norm)) %>%
  ungroup() -> sim_stats
target_stats$meansim = sim_stats$meansim[match(target_stats$gene, sim_stats$gene)]
target_stats$meansim[target_stats$n_launched_indic==1] = 1.00

write(paste("For all time, of ",nrow(target_stats)," launched targets, the ",
            sum(target_stats$n_launched_indic >= 10),
            " with ≥10 launched indications account for ",
            sum(target_stats$n_launched_indic[target_stats$n_launched_indic >= 10]),
            " of ",sum(target_stats$n_launched_indic),
            " launched indications (",
            percent(sum(target_stats$n_launched_indic[target_stats$n_launched_indic >= 10])/sum(target_stats$n_launched_indic)),
            ")",'\n',sep=''),text_stats_path,append=T)


# supplementary version with only year 2000+ programs
pp %>%
  filter(hcat=='Launched' & gene != '' & indication_mesh_id != '') %>%
  distinct(gene, indication_mesh_id) -> ti_launched
ti_launched %>%
  select(gene, meshcode_a = indication_mesh_id) -> a1
ti_launched %>%
  select(gene, meshcode_b = indication_mesh_id) -> a2
a1 %>%
  inner_join(a2, by=c('gene'='gene')) %>%
  filter(meshcode_a < meshcode_b) -> ti_simpairs
ti_simpairs %>%
  left_join(sim, by=c('meshcode_a' = 'meshcode_a', 'meshcode_b' = 'meshcode_b')) -> ti_sim
ti_sim$comb_norm[is.na(ti_sim$comb_norm)] = 0 # in current dataset this never arises - occurs only with corrupted mesh codes, which we fixed
ti_launched %>%
  group_by(gene) %>%
  summarize(.groups='keep', 
            n_launched_indic_2000 = length(unique(indication_mesh_id))) %>%
  ungroup() -> temp_target_stats
ti_sim %>%
  group_by(gene) %>%
  summarize(.groups='keep', 
            meansim_2000=mean(comb_norm)) %>%
  ungroup() -> temp_sim_stats
target_stats$n_launched_indic_2000 = temp_target_stats$n_launched_indic_2000[match(target_stats$gene, temp_target_stats$gene)]
target_stats$meansim_2000 = temp_sim_stats$meansim_2000[match(target_stats$gene, temp_sim_stats$gene)]
target_stats$meansim_2000[target_stats$n_launched_indic_2000==1] = 1.00

write(paste("Since 2000 only, of ",sum(target_stats$n_launched_indic_2000 > 0, na.rm=T)," launched targets, the ",
            sum(target_stats$n_launched_indic_2000 >= 10, na.rm=T),
            " with ≥10 launched indications account for ",
            sum(target_stats$n_launched_indic_2000[target_stats$n_launched_indic_2000 >= 10], na.rm=T),
            " of ",sum(target_stats$n_launched_indic_2000, na.rm=T),
            " launched indications (",
            percent(sum(target_stats$n_launched_indic_2000[target_stats$n_launched_indic_2000 >= 10], na.rm=T)/sum(target_stats$n_launched_indic_2000, na.rm=T)),
            ")",'\n',sep=''),text_stats_path,append=T)

hist_ti$ipert = target_stats$n_launched_indic[match(hist_ti$gene, target_stats$gene)]
hist_ti$meansim = target_stats$meansim[match(hist_ti$gene, target_stats$gene)]

hist_ti$ipert_2000 = target_stats$n_launched_indic_2000[match(hist_ti$gene, target_stats$gene)]
hist_ti$meansim_2000 = target_stats$meansim_2000[match(hist_ti$gene, target_stats$gene)]

areas$meansim = as.numeric(NA)
areas$meanipert = as.numeric(NA)

areas$meansim_2000 = as.numeric(NA)
areas$meanipert_2000 = as.numeric(NA)

areas$n_launched_ti = as.numeric(NA)
areas$n_gensup_ti = as.numeric(NA)

for (i in 1:nrow(areas)) {
  subset_by_area(hist_ti, areas$topl[i], areas$filter[i]) %>%
    select(ti_uid, gene, indication_mesh_id, ccat) %>%
    left_join(target_stats, by='gene') %>% 
    filter(ccat=='Launched') -> launched_area
  
  areas$meansim[i] = mean(launched_area$meansim)
  areas$meanipert[i] = mean(launched_area$n_launched_indic)
  areas$meansim_2000[i] = mean(launched_area$meansim_2000)
  areas$meanipert_2000[i] = mean(launched_area$n_launched_indic_2000)
  
  hist_ti_all_area = subset_by_area(hist_ti_all, areas$topl[i], 'only')
  areas$n_launched_ti[i] = length(unique(hist_ti_all_area$ti_uid[hist_ti_all_area$cat=='Launched']))
  areas$n_gensup_ti[i] = length(unique(hist_ti_all_area$ti_uid[hist_ti_all_area$target_status=='genetically supported target']))
}

nli_sim_spearman = suppressWarnings(cor.test(target_stats$n_launched_indic, target_stats$meansim, method='spearman'))
nli_sim_spearman_p = nli_sim_spearman$p.value
nli_sim_spearman_rho = nli_sim_spearman$estimate

write(paste("Spearman's correlation n_launched_indic vs. meansim across targets: rho = ",
            formatC(nli_sim_spearman_rho, format='fg', digits=3),
            ', P = ',formatC(nli_sim_spearman_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)


sim8 = sim[sim$comb_norm >= 0.8,]
assoc %>% filter(source != 'OTG' | l2g_share >= 0.5) -> assoc_l2g5
assoc_l2g5$associated_mesh_term = assoc_l2g5$mesh_term
indic %>% filter(genetic_insight != 'none') -> indic_insight
sim8 %>%
  inner_join(assoc_l2g5, by=c('meshcode_a'='mesh_id')) %>%
  inner_join(indic_insight, by=c('meshcode_b'='indication_mesh_id')) %>%
  group_by(gene, assoc_mesh_id=meshcode_a, indication_mesh_id=meshcode_b, source) %>%
  summarize(.groups='keep',
            assoc_mesh_term = min(associated_mesh_term),
            indication_mesh_term = min(indication_mesh_term),
            similarity = max(comb_norm),
            assoc_info = max(extra_info)) %>%
  ungroup() %>%
  mutate(ti_uid = paste0(gene,'-',indication_mesh_id)) -> all_possible_gensup_ti_assocs_all

read_tsv('data/areas.tsv', col_types=cols()) %>%
  select(topl, area, color) -> areas_temp

all_possible_gensup_ti_assocs_all %>%
  inner_join(indic_topl_match, by='indication_mesh_id') %>%
  inner_join(areas_temp, by='topl') %>%
  group_by(area) %>%
  summarize(.groups='keep', n_supp_ti = length(unique(ti_uid))) %>%
  ungroup() -> area_possible_supp_ti

areas$poss_supp_gi = area_possible_supp_ti$n_supp_ti[match(areas$area, area_possible_supp_ti$area)]





areas %>%
  select(topl, area, color, filter, rs_mean, rs_l95, rs_u95,
         nosup_clinical, nosup_launched, gensup_clinical, gensup_launched,
         p_s_mean, p_s_l95, p_s_u95, p_g_mean, p_g_l95, p_g_u95, poss_supp_gi,
         meansim, meanipert, meansim_2000, meanipert_2000, n_launched_ti, n_gensup_ti) -> area_output
write_supp_table(area_output, "Properties of therapy areas.")



ti_stats$cumnli = cumsum(ti_stats$n_launched_indic)
nti_tot = sum(ti_stats$n_launched_indic)
nti_tot_n = nrow(ti_stats)
nti_10plus = sum(ti_stats$n_launched_indic[ti_stats$n_launched_indic >= 10])
nti_10plus_n = sum(ti_stats$n_launched_indic >= 10)

hist_ti$gensup = hist_ti$target_status == 'genetically supported target'
logit_n_indic = glm(gensup ~ ipert, data=subset(hist_ti, cat=='Launched'), family='binomial')
# summary(logit_n_indic)
indic_per_target_beta = summary(logit_n_indic)$coefficients['ipert','Estimate']
indic_per_target_p    = summary(logit_n_indic)$coefficients['ipert','Pr(>|z|)']

write(paste("Logit model gensup ~ ipert: beta = ",
            formatC(indic_per_target_beta, format='fg', digits=3),
            ', P = ',formatC(indic_per_target_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)

logit_n_indic_no1 = glm(gensup ~ ipert, data=subset(hist_ti, cat=='Launched' & ipert > 1), family='binomial')
# summary(logit_n_indic)
indic_per_target_no1_beta = summary(logit_n_indic_no1)$coefficients['ipert','Estimate']
indic_per_target_no1_p    = summary(logit_n_indic_no1)$coefficients['ipert','Pr(>|z|)']

write(paste("Logit model gensup ~ ipert without 1-indication targets: beta = ",
            formatC(indic_per_target_no1_beta, format='fg', digits=3),
            ', P = ',formatC(indic_per_target_no1_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)

logit_meansim = glm(gensup ~ meansim, data=subset(hist_ti, cat=='Launched'), family='binomial')
# summary(logit_meansim)
meansim_beta = summary(logit_meansim)$coefficients['meansim','Estimate']
meansim_p = summary(logit_meansim)$coefficients['meansim','Pr(>|z|)']

write(paste("Logit model gensup ~ meansim: beta = ",
            formatC(meansim_beta, format='fg', digits=3),
            ', P = ',formatC(meansim_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)

logit_meansim_no1 = glm(gensup ~ meansim, data=subset(hist_ti, cat=='Launched' & ipert > 1), family='binomial')
# summary(logit_meansim)
meansim_no1_beta = summary(logit_meansim_no1)$coefficients['meansim','Estimate']
meansim_no1_p = summary(logit_meansim_no1)$coefficients['meansim','Pr(>|z|)']

write(paste("Logit model gensup ~ meansim without 1-indication targets: beta = ",
            formatC(meansim_no1_beta, format='fg', digits=3),
            ', P = ',formatC(meansim_no1_p, format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)

resx=300
png(paste0(output_path,'/figure-2.png'),width=6.5*resx,height=8.5*resx,res=resx)

#layout_matrix =  matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,6,
#                          7,7,7,8,8,8,8,9,9,9,10,10,10), nrow=2, byrow=T)

layout_matrix =  matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,6,
                          7,7,7,7,7,8,8,8,8,9,9,9,9,
                          10,10,10,10,10,11,11,11,11,12,12,12,12), nrow=3, byrow=T)

layout(layout_matrix, heights=c(1,.65,.65))

line_color = '#B9B9B9'
xlims = c(0,5)
ylims = range(areas_all$y, na.rm=T) + c(-0.5, 0.5)
par(mar=c(3,1,3,0))
plot(NA, NA, ylim=ylims, xlim=xlims, axes=F, ann=F, xaxs='i', yaxs='i')
areas_all %>% group_by(topl, area, color, y) %>% slice(1) %>% select(topl, area, color, y) -> rr_table_meta
mtext(side=4, at=rr_table_meta$y, line=0, adj=1, las=2, text=rr_table_meta$area, cex=.8, col=rr_table_meta$color)
panel = 1
par(mar=c(3,0.25,3,3.5))
transition_disp = list('Preclinical'='Pre-I','I'='I-II', 'II'='II-III', 'III'='III-Launch', 'I-Launch'='I-Launch')
for (transition in c('Preclinical','I','II','III','I-Launch')) {
  subs = subset(rr_table_areas, phase==transition & !is.na(y))
  if (transition == 'I-Launch') {
    xlims = c(0, 5)
  } else {
    xlims = c(0, 2.5)
  }
  plot(NA, NA, ylim=ylims, xlim=xlims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
  axis(side=1, lwd=0, lwd.ticks=1, at=0:100/10, labels=NA, tck=-0.015)
  axis(side=1, lwd=0, lwd.ticks=1, at=0:4, labels=NA, tck=-0.03)
  axis(side=1, lwd=0, at=0:5, line=-0.5)
  mtext(side=1, line=2, text='RS')
  abline(v=0:4, lwd=0.125, col=line_color)
  abline(v=1, lwd=.375, lty=3, col='black')
  axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0)
  mtext(side=4, at=subs$y, line=0.5, las=2, text=subs$fraction, cex=0.6, col=subs$color)
  segments(x0=subs$rs_l, x1=pmin(subs$rs_u,max(xlims)+1), y0=subs$y, lwd=2, col=subs$color)
  points(x=subs$rs_mean, y=subs$y, pch=19, col=subs$color)
  mtext(side=3, line=0.25, text=transition_disp[[transition]], adj=0, at=1, cex=0.9, font=2)
  
  par(xpd=T)
  segments(x0=-max(xlims), x1=max(xlims)*1.25, y0=max(areas_all$y)-0.5, lwd=.125)
  #segments(x0=-max(xlims), x1=max(xlims)*1.25, y0=0, lwd=.125)
  par(xpd=F)
  
  mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
  panel = panel + 1
}

par(mar=c(4,4,3,1))
indic %>% filter(genetic_insight != 'none') -> indic_insight
sim %>%
  filter(comb_norm >= 0.8) %>%
  filter(meshcode_a %in% indic_insight$indication_mesh_id) -> sim8

assoc %>%
  filter(source == 'OTG' & l2g_share >= 0.5) %>%
  group_by(gene, mesh_id) %>%
  summarize(.groups='keep', min_year = min(year)) %>%
  ungroup() -> gty

# take only entries where mesh code A is a genetic insight indication, 
# i.e. ever developed in pharmaprojects and has genetic insight
# (this also excludes diagnostic indications)
gty %>%
  inner_join(sim8, by=c('mesh_id'='meshcode_b')) %>%
  group_by(gene, meshcode_a) %>%
  summarize(.groups='keep', min_year = min(min_year)) %>%
  ungroup() -> tiy

gty %>%
  group_by(min_year) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() %>%
  arrange(min_year) %>%
  mutate(cumn = cumsum(n)) -> gtyc

tiy %>%
  group_by(min_year) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() %>%
  arrange(min_year) %>%
  mutate(cumn = cumsum(n)) -> tiyc

gtyc %>%
  inner_join(tiyc, by='min_year', suffix=c('_gt','_ti')) -> gtyc_tiyc

write_supp_table(gtyc_tiyc, 'Cumulative OTG G-T associations and supported T-I pairs by year.')

read_tsv('data/areas.tsv', col_types=cols()) %>%
  select(topl, area, color) -> areas_temp

tiy %>%
  inner_join(indic_topl_match, by=c('meshcode_a'='indication_mesh_id')) %>%
  group_by(min_year, topl) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() %>%
  arrange(min_year, topl) %>%
  group_by(topl) %>%
  mutate(cumn = cumsum(n)) %>%
  ungroup() %>%
  inner_join(areas_temp, by='topl') -> tiya

tiya %>%
  group_by(topl, area, color) %>%
  slice_max(min_year) %>%
  ungroup() %>%
  arrange(desc(cumn)) -> tiya_leg

par(mar=c(3,4,2,6))
xlims = c(2007, 2022)
xats = 2007:2022
xbigs = c(2010, 2015, 2020)
ylims = c(0, 27000)
yats = 0:10*1e3
ybigs = c(0, 5000, 10000, 15000, 20000, 25000)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, line=-0.5, lwd=0, labels=xbigs)
mtext(side=1, line=2.0, text='Year')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
axis(side=2, at=yats, tck=-0.015, labels=NA)
axis(side=2, at=ybigs, tck=-0.03, labels=NA)
axis(side=2, at=ybigs, line=-0.5, lwd=0, labels=paste0(ybigs/1000,'K'), las=2)
mtext(side=2, line=2, text='Cumulative G-I pairs')
for (this_area in unique(tiya$area)) {
  subs = tiya[tiya$area==this_area,]
  points(subs$min_year, subs$cumn, col=subs$color, lwd=1.5, type='l')
}
par(xpd=T)
legend(x=max(xlims), y=max(ylims)*1.1, tiya_leg$area, col=tiya_leg$color, text.col=tiya_leg$color, lwd=1.5, bty='n', cex=0.7)
par(xpd=F)
mtext(letters[panel], side=3, cex=2, adj = 0.0, line = 0.5)
panel = panel + 1

write_supp_table(tiya, 'Cumulative OTG supported G-I pairs by year and therapy area.')

cex_factor = 70
par(mar=c(3,2.5,2,1))
xlims = c(0,30000)
plot(NA, NA, xlim=xlims, ylim=c(0,4.5), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, lwd=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, lwd=0, lwd.ticks=1, at=0:30*1e3, labels=NA, tck=-0.015)
axis(side=1, lwd=0, lwd.ticks=1, at=0:6*5e3, labels=NA, tck=-0.03)
axis(side=1, lwd=0, line=-0.5, at=0:6*5e3, labels=c('0',paste0(1:6*5,'K')))
mtext(side=1, line=1.75, text="Possible supported G-I")
mtext(side=2, line=1.5, text='RS')
axis(side=2, lwd=1, lwd.ticks=1, at=0:10/2, labels=NA, tck=-0.015)
axis(side=2, lwd=0, lwd.ticks=1, at=0:5, labels=NA, tck=-0.03)
axis(side=2, lwd=0, at=0:5, line=-0.5, las=2)
abline(h=0:5, lwd=0.125, col=line_color)
abline(h=1, lwd=0.5, col='black')
par(xpd=T)
points(x=areas$poss_supp_gi, y=areas$rs_mean, col=areas$color, pch=19, cex=areas$n_gensup_ti/cex_factor)
par(xpd=F)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1
allsupp_rs_weightedbygensup = wtd.cor(areas$poss_supp_gi, areas$rs_mean, weight=areas$n_gensup_ti)
allsupp_rs_weightedbygensup_rho = round(allsupp_rs_weightedbygensup['Y','correlation'],2)
allsupp_rs_weightedbygensup_p = round(allsupp_rs_weightedbygensup['Y','p.value'],3)

write(paste('Weighted Pearson correlation across areas, poss_supp_gi vs. rs: rho= ',formatC(allsupp_rs_weightedbygensup_rho,digits=2,format='fg'),', P = ',formatC(allsupp_rs_weightedbygensup_p,digits=2,format='e'),'\n',sep=''),text_stats_path,append=T)



par(mar=c(3,4,2,1))
plot(NA, NA, xlim=c(0,55), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, lwd=1, lwd.ticks=1, at=0:100, labels=NA, tck=-0.015)
axis(side=1, lwd=0, lwd.ticks=1, at=0:10*10, labels=NA, tck=-0.03)
axis(side=1, lwd=0, at=0:10*10, line=-0.5)
mtext(side=1, line=1.75, text='Approved indications')
mtext(side=2, line=2.5, text='Mean similarity')
axis(side=2, lwd=1, lwd.ticks=1, at=0:100/100, labels=NA, tck=-0.015)
axis(side=2, lwd=0, lwd.ticks=1, at=0:10/10, labels=NA, tck=-0.03)
axis(side=2, lwd=0, at=0:10/10, line=-0.25, las=2)
set.seed(1)
points(target_stats$n_launched_indic, target_stats$meansim, pch=20)
callouts = tibble(gene = c('NR3C1','PTGS2','IFNAR2','TOP2A','OPRM1','PDCD1','CHRM3','SLC6A4'),
                      pos=c(4, 4, 3, 4, 4, 3, 1, 4))
callouts$x = target_stats$n_launched_indic[match(callouts$gene, target_stats$gene)]
callouts$y = target_stats$meansim[match(callouts$gene, target_stats$gene)]
callouts$x_orig = target_stats$n_launched_indic[match(callouts$gene, target_stats$gene)]
callouts$y_orig = target_stats$meansim[match(callouts$gene, target_stats$gene)]
callouts$x = callouts$x_orig
callouts$y = callouts$y_orig
callouts$x[callouts$gene %in% c('SLC6A4')] = callouts$x_orig[callouts$gene %in% c('SLC6A4')] + 8
callouts$y[callouts$gene %in% c('SLC6A4')] = callouts$y_orig[callouts$gene %in% c('SLC6A4')] - 0.13
callouts$x[callouts$gene %in% c('CHRM3')] = callouts$x_orig[callouts$gene %in% c('CHRM3')] - 2
callouts$y[callouts$gene %in% c('CHRM3')] = callouts$y_orig[callouts$gene %in% c('CHRM3')] - 0.05
callouts$x[callouts$gene %in% c('PTGS2')] = callouts$x_orig[callouts$gene %in% c('PTGS2')] + 4
callouts$y[callouts$gene %in% c('PTGS2')] = callouts$y_orig[callouts$gene %in% c('PTGS2')] - 0.04
callouts$x[callouts$gene %in% c('OPRM1')] = callouts$x_orig[callouts$gene %in% c('OPRM1')] + 6
callouts$y[callouts$gene %in% c('OPRM1')] = callouts$y_orig[callouts$gene %in% c('OPRM1')] - 0.15
callouts$x[callouts$gene %in% c('TOP2A')] = callouts$x_orig[callouts$gene %in% c('TOP2A')] + 5
callouts$y[callouts$gene %in% c('TOP2A')] = callouts$y_orig[callouts$gene %in% c('TOP2A')] + 0.1
callouts$x[callouts$gene %in% c('PDCD1')] = callouts$x_orig[callouts$gene %in% c('PDCD1')] + 0
callouts$y[callouts$gene %in% c('PDCD1')] = callouts$y_orig[callouts$gene %in% c('PDCD1')] + 0.1
par(xpd=T)
segments(x0=callouts$x_orig, x1=callouts$x, y0=callouts$y_orig, y1=callouts$y, lwd=0.75)
text(x=callouts$x, y=callouts$y, pos=callouts$pos, font=3, labels=callouts$gene, cex=0.8)
par(xpd=F)
indic_sim_spearman = suppressWarnings(cor.test(target_stats$n_launched_indic, target_stats$meansim, method='spearman'))
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(target_stats, "Count and mean similarity of approved indications per target.")

write(paste('Spearman correlation n_launched_indic vs. mean_sim: rho= ',formatC(indic_sim_spearman$estimate,digits=2,format='fg'),', P = ',formatC(indic_sim_spearman$p.value,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)

indic_per_quantiles = quantile(hist_ti$ipert[hist_ti$cat=='Launched'], 0:5/5, na.rm=T)
indic_per_forest = data.frame(y=5:1, 
                              min=c(1,floor(indic_per_quantiles[c('20%','40%','60%','80%')])),
                              max=c(floor(indic_per_quantiles[c('20%','40%','60%','80%')])-1,round(indic_per_quantiles['100%']))) # use >= min and <= max
# note - there are no non-integer values in indic_per_target
# which(hist_ti$indic_per_target[hist_ti$cat=='Launched'] != round(hist_ti$indic_per_target[hist_ti$cat=='Launched'],0))
# integer(0)
# but if you do quantile(hist_ti$indic_per_target[hist_ti$cat=='Launched'], .4) you get 7.8
# because quantile() uses type 7 algorithm by default which gives an interpolated average if the numbers on either side are different
for (i in 1:nrow(indic_per_forest)) {
  all_launched = length(unique(hist_ti$ti_uid[hist_ti$cat=='Launched' & hist_ti$ipert >= indic_per_forest$min[i] & hist_ti$ipert <= indic_per_forest$max[i]]))
  gensup_launched = length(unique(hist_ti$ti_uid[hist_ti$cat=='Launched' & hist_ti$ipert >= indic_per_forest$min[i] & hist_ti$ipert <= indic_per_forest$max[i] & hist_ti$gensup]))
  bconf_obj = binom.confint(x=gensup_launched, n=all_launched, method='wilson')
  indic_per_forest[i,'label'] = ifelse(indic_per_forest$min[i]==indic_per_forest$max[i], indic_per_forest$min[i], paste0(indic_per_forest$min[i],'-',indic_per_forest$max[i]))
  indic_per_forest$numerator[i] = gensup_launched
  indic_per_forest$denominator[i] = all_launched
  indic_per_forest[i,c('mean','l95','u95')] = bconf_obj[,c('mean','lower','upper')]
}

meansim_forest = tibble(min=seq(0,0.99,0.2),
                        max=c(seq(0.2,0.99,0.2),1.01)) %>% # use >= min and < max, hence top bound is 1.01 not 1.0
  mutate(y = max(row_number()) - row_number() + 1) %>%
  mutate(label = paste0(formatC(min,format='f',digits=1),'-',substr(as.character(max),1,3))) %>%
  mutate(numerator=as.numeric(NA), 
         denominator=as.numeric(NA),
         mean=as.numeric(NA),
         l95=as.numeric(NA),
         u95=as.numeric(NA))
for (i in 1:nrow(meansim_forest)) {
  all_launched = length(unique(hist_ti$ti_uid[hist_ti$cat=='Launched' & hist_ti$meansim >= meansim_forest$min[i] & hist_ti$meansim < meansim_forest$max[i]]))
  gensup_launched = length(unique(hist_ti$ti_uid[hist_ti$cat=='Launched' & hist_ti$meansim >= meansim_forest$min[i] & hist_ti$meansim < meansim_forest$max[i] & hist_ti$gensup]))
  bconf_obj = binom.confint(x=gensup_launched, n=all_launched, method='wilson')
  meansim_forest$numerator[i] = gensup_launched
  meansim_forest$denominator[i] = all_launched
  meansim_forest[i,c('mean','l95','u95')] = bconf_obj[,c('mean','lower','upper')]
}

master_forest_2 = rbind(indic_per_forest, meansim_forest)
master_forest_2$y = c(11:7,5:1)


plot_forest(master_forest_2, xlims=c(0,.35), xstyle='percent', mar=c(4,8,2.5,6))
tranche_line = 4.0
axis(side=2, at=c(6.75,11.25), tck=0.025, labels=NA, line=tranche_line)
axis(side=2, at=c(0.75,5.25), tck=0.025, labels=NA, line=tranche_line)
abline(h=6)
mtext(side=2, at=c(3.5, 9.5), text=c('Mean\nsimilarity', 'Indications/\ntarget'), line=tranche_line + 0.25, cex=0.8)
mtext(side=1, line=2, text='P(G)')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(master_forest_2, "Probability of genetic support by count and similarity of indications per target.")

cex_factor = 70

par(mar=c(4,3,2.5,1))
plot(NA, NA, xlim=c(0.3,0.65), ylim=c(0,4.5), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, lwd=1, at=c(0,1), lwd.ticks=0, labels=NA)
axis(side=1, lwd=0, lwd.ticks=1, at=0:20/20, labels=NA, tck=-0.015)
axis(side=1, lwd=0, lwd.ticks=1, at=0:10/10, labels=NA, tck=-0.03)
axis(side=1, lwd=0, at=0:10/10, line=-0.5)
mtext(side=1, line=1.5, text="Mean similarity")
mtext(side=2, line=1.75, text='RS')
axis(side=2, lwd=1, lwd.ticks=1, at=0:10/2, labels=NA, tck=-0.015)
axis(side=2, lwd=0, lwd.ticks=1, at=0:5, labels=NA, tck=-0.03)
axis(side=2, lwd=0, at=0:5, line=-0.5, las=2)
abline(h=0:5, lwd=0.125, col=line_color)
abline(h=1, lwd=0.5, col='black')
par(xpd=T)
points(x=areas$meansim, y=areas$rs_mean, col=areas$color, pch=19, cex=areas$n_gensup_ti/cex_factor)
par(xpd=F)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1
sim_rr_weightedbygensup = wtd.cor(areas$meansim, areas$rs_mean, weight=areas$n_gensup_ti)
sim_rr_weightedbygensup_rho = round(sim_rr_weightedbygensup['Y','correlation'],2)
sim_rr_weightedbygensup_p = round(sim_rr_weightedbygensup['Y','p.value'],3)

write(paste('Weighted Pearson correlation across areas, mean_sim vs. rs: rho= ',formatC(sim_rr_weightedbygensup_rho,digits=2,format='fg'),', P = ',formatC(sim_rr_weightedbygensup_p,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)

par(mar=c(4,3,2.5,1))
xlims = c(0,18)
plot(NA, NA, xlim=xlims, ylim=c(0,4.5), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, lwd=0, lwd.ticks=1, at=0:100, labels=NA, tck=-0.015)
axis(side=1, lwd=0, lwd.ticks=1, at=0:10*10, labels=NA, tck=-0.03)
axis(side=1, lwd=0, at=0:10*10, line=-0.5)
mtext(side=1, line=1.5, text='Indications/target')
mtext(side=2, line=1.75, text='RS')
axis(side=2, lwd=1, lwd.ticks=1, at=0:10/2, labels=NA, tck=-0.015)
axis(side=2, lwd=0, lwd.ticks=1, at=0:5, labels=NA, tck=-0.03)
axis(side=2, lwd=0, at=0:5, line=-0.5, las=2)
abline(h=0:5, lwd=0.125, col=line_color)
abline(h=1, lwd=0.5, col='black')
par(xpd=T)
points(x=areas$meanipert, y=areas$rs_mean, col=areas$color, pch=19, cex=areas$n_gensup_ti/cex_factor)
par(xpd=F)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1
ipert_rr_weightedbygensup = wtd.cor(areas$meanipert, areas$rs_mean, weight=areas$n_gensup_ti)
ipert_rr_weightedbygensup_rho = round(ipert_rr_weightedbygensup['Y','correlation'],2)
ipert_rr_weightedbygensup_p = round(ipert_rr_weightedbygensup['Y','p.value'],3)

write(paste('Weighted Pearson correlation across areas, indic per target vs. rs: rho= ',formatC(ipert_rr_weightedbygensup_rho,digits=2,format='fg'),', P = ',formatC(ipert_rr_weightedbygensup_p,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)

unecessary_message = dev.off()











resx=300
png('display_items/figure-s5.png',width=6.5*resx,height=8.0*resx,res=resx)
panel = 1
par(mfrow = c(6,3))
for (i in 1:nrow(areas_all)) {
  
  combined_ti_area = subset_by_area(combined_ti, topl=areas_all$topl[i], filter=areas_all$filter[i])
  area_forest = advancement_forest(combined_ti_area)
  left_margin = case_when(i %% 3 == 1 ~ 6,
                          TRUE ~ 0)
  bottom_margin = case_when(i >= 16 ~ 3,
                            TRUE ~ 0)
  margins = c(bottom_margin, left_margin, 2, 0)
  
  plot_forest(area_forest, xlims=c(0,.60), xlab='', col=areas_all$color[i], title='', mar=margins)
  mtext(side=3, adj=0, at=-0.3, line=0.25, col='#000000', text=areas_all$area[i], cex=.85)
  # mtext(letters[panel], side=3, cex=1.4, at=-0.5, line = 0.25)
  # panel = panel + 1
  
  area_forest$area = areas_all$area[i]
  if (i == 1) {
    pg_out = area_forest
  } else {
    pg_out = rbind(pg_out, area_forest)
  }
}
unnecessary_message = dev.off()

pg_out %>%
  select(area, phase=label, phasenum=num, supported=numerator, total=denominator, pg_mean=mean, pg_l95 = l95, pg_u95 = u95) -> pg_out

write_supp_table(pg_out, "Proportion of target-indication pairs with genetic support by phase and therapy area, combined mode (both historical and active programs).")









cat(file=stderr(), 'done.\nCreating Figure S6...')

resx=300
png(paste0(output_path,'/figure-s6.png'),width=6.5*resx,height=9*resx,res=resx)
layout_matrix = matrix(c(1,1,2,2,2,3,3,3,9,
                         4,4,4,5,5,5,6,6,6,
                         7,7,7,7,7,8,8,8,8,
                         7,7,7,7,7,8,8,8,8), nrow=4,byrow=T)
layout(layout_matrix, heights=c(1.5,1,.75,.75))
panel = 1

ylims = range(areas$y, na.rm=T) + c(-0.5, 0.5)
par(mar=c(4,0,3,0))
plot(NA, NA, ylim=ylims, xlim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')
mtext(side=4,adj=1,las=2,at=areas$y,text=areas$area,col=areas$color)

xlims = c(0, 0.20)
par(mar=c(4,1,3,8))
plot(NA, NA, ylim=ylims, xlim=xlims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, lwd=0, lwd.ticks=1, at=0:100/100, labels=NA, tck=-0.015)
axis(side=1, lwd=0, lwd.ticks=1, at=0:10/10, labels=NA, tck=-0.03)
axis(side=1, lwd=0, at=0:10/10, labels=percent(0:10/10), line=-0.5)
mtext(side=1, line=2, text='P(S)')
abline(v=0:4, lwd=0.125, col=line_color)
abline(v=1, lwd=1, col='black')
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0)
mtext(side=4, at=areas$y, line=0.5, las=2, text=paste0(formatC(areas$gensup_launched+areas$nosup_launched, format='d', big.mark=','),'/',formatC(areas$nosup_launched+areas$gensup_launched+areas$nosup_clinical+areas$gensup_clinical, format='d', big.mark=',')), cex=1, col=areas$color)
segments(x0=areas$p_s_l95, x1=pmin(areas$p_s_u95,max(xlims)+1), y0=areas$y, lwd=1.5, col=areas$color)
points(x=areas$p_s_mean, y=areas$y, pch=19, col=areas$color)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

ylims = range(areas$y, na.rm=T) + c(-0.5, 0.5)
xlims = c(0, 0.30)
par(mar=c(4,1,3,8))
plot(NA, NA, ylim=ylims, xlim=xlims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, lwd=0, lwd.ticks=1, at=0:100/100, labels=NA, tck=-0.015)
axis(side=1, lwd=0, lwd.ticks=1, at=0:10/10, labels=NA, tck=-0.03)
axis(side=1, lwd=0, at=0:10/10, labels=percent(0:10/10), line=-0.5)
mtext(side=1, line=2, text='P(G)')
abline(v=0:4, lwd=0.125, col=line_color)
abline(v=1, lwd=1, col='black')
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0)
mtext(side=4, at=areas$y, line=0.5, las=2, text=paste0(formatC(areas$gensup_launched+areas$gensup_clinical, format='d', big.mark=','),'/',formatC(areas$nosup_launched+areas$gensup_launched+areas$nosup_clinical+areas$gensup_clinical, format='d', big.mark=',')), cex=1, col=areas$color)
segments(x0=areas$p_g_l95, x1=pmin(areas$p_g_u95,max(xlims)+1), y0=areas$y, lwd=1.5, col=areas$color)
points(x=areas$p_g_mean, y=areas$y, pch=19, col=areas$color)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

ylims = c(0.0, 0.20)
xlims = c(0.0, 0.30)
par(mar=c(4,5,3,1))
plot(NA, NA, ylim=ylims, xlim=xlims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, lwd=0, lwd.ticks=1, at=0:100/100, labels=NA, tck=-0.015)
axis(side=1, lwd=0, lwd.ticks=1, at=0:10/10, labels=NA, tck=-0.03)
axis(side=1, lwd=0, at=0:10/10, labels=percent(0:10/10), line=-0.5)
mtext(side=1, line=2, text='P(G)')
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=2, lwd=0, lwd.ticks=1, at=0:100/100, labels=NA, tck=-0.015)
axis(side=2, lwd=0, lwd.ticks=1, at=0:10/10, labels=NA, tck=-0.03)
axis(side=2, lwd=0, at=0:10/10, labels=percent(0:10/10), line=-0.5, las=2)
mtext(side=2, line=2.5, text='P(S)')
segments(x0=areas$p_g_l95, x1=pmin(areas$p_g_u95,max(xlims)+1), y0=areas$p_s_mean, lwd=1.5, col=areas$color)
segments(x0=areas$p_g_mean, y0=areas$p_s_l95, y1=areas$p_s_u95, lwd=1.5, col=areas$color)
points(x=areas$p_g_mean, y=areas$p_s_mean, pch=19, col=areas$color, cex=1.5)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

ylims = c(0,5)
xlims = c(0.0, 0.20)
par(mar=c(4,5,3,1))
plot(NA, NA, ylim=ylims, xlim=xlims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, lwd=0, lwd.ticks=1, at=0:100/100, labels=NA, tck=-0.015)
axis(side=1, lwd=0, lwd.ticks=1, at=0:10/10, labels=NA, tck=-0.03)
axis(side=1, lwd=0, at=0:10/10, labels=percent(0:10/10), line=-0.5)
mtext(side=1, line=2, text='P(S)')
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=2, lwd=0, lwd.ticks=1, at=0:100/10, labels=NA, tck=-0.015)
axis(side=2, lwd=0, lwd.ticks=1, at=0:10, labels=NA, tck=-0.03)
axis(side=2, lwd=0, at=0:10, labels=0:10, line=-0.5, las=2)
abline(h=1,lwd=0.5)
mtext(side=2, line=2.5, text='RS')
segments(x0=areas$p_s_l95, x1=pmin(areas$p_s_u95,max(xlims)+1), y0=areas$rs_mean, lwd=1.5, col=areas$color)
segments(x0=areas$p_s_mean, y0=areas$rs_l95, y1=areas$rs_u95, lwd=1.5, col=areas$color)
points(x=areas$p_s_mean, y=areas$rs_mean, pch=19, col=areas$color, cex=1.5)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

ylims = c(0,5)
xlims = c(0, 0.30)
par(mar=c(4,5,3,1))
plot(NA, NA, ylim=ylims, xlim=xlims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, lwd=0, lwd.ticks=1, at=0:100/100, labels=NA, tck=-0.015)
axis(side=1, lwd=0, lwd.ticks=1, at=0:10/10, labels=NA, tck=-0.03)
axis(side=1, lwd=0, at=0:10/10, labels=percent(0:10/10), line=-0.5)
mtext(side=1, line=2, text='P(G)')
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=2, lwd=0, lwd.ticks=1, at=0:100/10, labels=NA, tck=-0.015)
axis(side=2, lwd=0, lwd.ticks=1, at=0:10, labels=NA, tck=-0.03)
axis(side=2, lwd=0, at=0:10, labels=0:10, line=-0.5, las=2)
abline(h=1,lwd=0.5)
mtext(side=2, line=2.5, text='RS')
segments(x0=areas$p_g_l95, x1=pmin(areas$p_g_u95,max(xlims)+1), y0=areas$rs_mean, lwd=1.5, col=areas$color)
segments(x0=areas$p_g_mean, y0=areas$rs_l95, y1=areas$rs_u95, lwd=1.5, col=areas$color)
points(x=areas$p_g_mean, y=areas$rs_mean, pch=19, col=areas$color, cex=1.5)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1
cor_ps_pg_wtd = wtd.cor(areas$p_g_mean, areas$p_s_mean, w=areas$n_gensup_ti)
cor_ps_rs_wtd = wtd.cor(areas$p_s_mean, areas$rs_mean, w=areas$n_gensup_ti)
cor_pg_rs_wtd = wtd.cor(areas$p_g_mean, areas$rs_mean, w=areas$n_gensup_ti)
cor_pg_rs_unwtd = cor.test(areas$p_g_mean, areas$rs_mean, method='pearson')
cor_ps_pg_unwtd = cor.test(areas$p_s_mean, areas$p_g_mean, method='pearson')
cor_ps_rs_unwtd = cor.test(areas$p_s_mean, areas$rs_mean, method='pearson')
model_ps_pg_rs = lm(rs_mean ~ p_g_mean + p_s_mean, data=areas)

write(paste("Weighted Pearson's correlation P(G) vs P(S) across areas: rho = ",
            formatC(cor_ps_pg_wtd['Y','correlation'], format='fg', digits=3),
            ', P = ',formatC(cor_ps_pg_wtd['Y','p.value'], format='e', digits=1),
            '\n',sep=''),text_stats_path,append=T)

write(paste("Weighted Pearson's correlation P(S) vs RS across areas: rho = ",
            formatC(cor_ps_rs_wtd['Y','correlation'], format='fg', digits=3),
            ', P = ',formatC(cor_ps_rs_wtd['Y','p.value'], format='e', digits=1),
      '\n',sep=''),text_stats_path,append=T)

write(paste("Weighted Pearson's correlation P(G) vs RS across areas: rho = ",
            formatC(cor_pg_rs_wtd['Y','correlation'], format='fg', digits=3),
            ', P = ',formatC(cor_pg_rs_wtd['Y','p.value'], format='e', digits=1),
      '\n',sep=''),text_stats_path,append=T)

tbl_s6 = read_tsv('data/nelson_2015_table_s6.tsv', col_types=cols())
colnames(tbl_s6) = gsub('[^a-z0-9_]','_',tolower(colnames(tbl_s6)))
tbl_s6$ti = paste0(tbl_s6$gene, '-', tbl_s6$msh_ind)
phase_mapping = data.frame(orig=unique(tbl_s6$latestphase))
phase_mapping$order = c('4 Approved','0 Preclinical','1 Phase I','2 Phase II','3 Phase III')
tbl_s6$phase = phase_mapping$order[match(tbl_s6$latestphase, phase_mapping$orig)]
all_cat = data.frame(category=sort(unique(tbl_s6$category)))
tbl_s6 %>%
  group_by(ti, gene, msh_ind, category) %>%
  summarize(.groups='keep', maxphase = max(phase)) -> gensup_ti
all_cat_phase = expand.grid(phase=phase_mapping$order, category=all_cat$category)
left_join(all_cat_phase, gensup_ti, by=c('phase'='maxphase', 'category'='category')) %>%
  group_by(phase, category) %>%
  summarize(.groups='keep', n_ti_gensup = length(unique(ti))) -> gensup_cat
# pasted into Excel and manually added in totals from Fig 3B & S4A
data2015 = read_tsv('data/nelson_2015_by_category.tsv', col_types=cols())
match_2015_2020 = read_tsv('data/nelson_2015_category_match.tsv', col_types=cols())
data2015 %>%
  mutate(n_ti_nosup = n_ti_total - n_ti_gensup) %>%
  group_by(category) %>%
  summarize(.groups='keep',
            nosup_clinical = sum(n_ti_nosup * (phase %in% c('1 Phase I','2 Phase II','3 Phase III'))),
            nosup_launched = sum(n_ti_nosup * (phase %in% c('4 Approved'))),
            gensup_clinical = sum(n_ti_gensup * (phase %in% c('1 Phase I','2 Phase II','3 Phase III'))),
            gensup_launched = sum(n_ti_gensup * (phase %in% c('4 Approved')))) -> smry2015
smry2015$rs_mean = as.numeric(NA)
smry2015$rs_l95 = as.numeric(NA)
smry2015$rs_u95 = as.numeric(NA)
for (i in 1:nrow(smry2015)) {
  ctable = matrix(as.integer(smry2015[i,c("nosup_clinical", "nosup_launched", "gensup_clinical", "gensup_launched")]), nrow=2, byrow=T, dimnames=list(gensup=c('0','1'),phase=c('clinical','launched')))
  rr_epitools_obj = suppressWarnings(riskratio(ctable, method='boot', replicates=50000))
  smry2015[i,c('rs_mean','rs_l95','rs_u95')] = as.list(rr_epitools_obj$measure['1',c('estimate','lower','upper')])
}
smry2015 = smry2015[with(smry2015, order(-rs_mean)),]
smry2015$y = rank(smry2015$rs_mean, ties.method = 'first')
smry2015$color = '#000000'
smry2015$area = smry2015$category
ylims = range(smry2015$y, na.rm=T) + c(-0.5, 0.5)
smry2015$rs_u95[smry2015$rs_u95==Inf] = max(ylims)
xlims = c(0, 8)
par(mar=c(4,10,3,4))
plot(NA, NA, ylim=ylims, xlim=xlims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, lwd=0, lwd.ticks=1, at=0:100/10, labels=NA, tck=-0.015, cex.axis=0.7)
axis(side=1, lwd=0, lwd.ticks=1, at=0:10, labels=NA, tck=-0.03, cex.axis=0.7)
axis(side=1, lwd=0, at=0:10, line=-0.5)
abline(v=0:8, lwd=0.125, col=line_color)
abline(v=1, lwd=1, col='black')
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0)
mtext(side=2, at=smry2015$y, line=0.5, las=2, cex=0.7, text=smry2015$area, col=smry2015$color)
mtext(side=4, at=smry2015$y, line=0.5, las=2, cex=0.7, text=paste0(smry2015$gensup_launched,'/',smry2015$gensup_launched + smry2015$gensup_clinical), col=smry2015$color)
segments(x0=smry2015$rs_l95, x1=smry2015$rs_u95, y0=smry2015$y, lwd=3, col=smry2015$color)
points(x=smry2015$rs_mean, y=smry2015$y, pch=19, col=smry2015$color)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1
smry2015$n_total = rowSums(smry2015[,c("nosup_clinical", "nosup_launched", "gensup_clinical", "gensup_launched")])
wmean2015 = round(weighted.mean(smry2015$rs_mean, w=smry2015$n_total),2)
wmean2022 = round(weighted.mean(areas$rs_mean, w=areas$n_hist_ti),2)
ctable = matrix(as.integer(colSums(smry2015[,c("nosup_clinical", "nosup_launched", "gensup_clinical", "gensup_launched")])), nrow=2, byrow=T, dimnames=list(gensup=c('0','1'),phase=c('clinical','launched')))
rr_epitools_obj = riskratio(ctable, method='boot', replicates=50000)
total2015 = round(rr_epitools_obj$measure['1',c('estimate','lower','upper')],2)
rr_obj = advancement_rr(hist_ti)
total2022 = round(rr_obj[rr_obj$phase=='I-Launch',c('rs_mean','rs_l','rs_u')],2)

write_supp_table(smry2015, "Summary of therapy area data from Nelson 2015")

write(paste("Overall RS from Nelson 2015: ",
            formatC(total2015['estimate'], format='fg', digits=3),
            ' vs. current study: ',formatC(total2022$rs_mean, format='fg', digits=3),
      '\n',sep=''),text_stats_path,append=T)
# 
# # cross-tab of 2013 area assignments vs. 2023 area assignments
# n15p = read_tsv('data/nelson_2015_supplementary_dataset_3.tsv', col_types=cols()) %>%
#   clean_names() %>%
#   filter(!(phase_latest %in% c('Discontinued','No Development Reported'))) %>%
#   mutate(ccat = gsub(' Clinical Trial','',phase_latest)) %>%
#   mutate(mesh = tolower(msh))
# n15p$mesh_id = vocab_match$id[match(n15p$msh, vocab_match$labeltext)]
# n15p$mesh_term = mesh_best_names$labeltext[match(n15p$mesh_id, mesh_best_names$id)]

# TO DO: table s6 is only the genetically supported ones.
# instead need to switch to supp dataset 3 which is all of PP in 2013
# then use supp table 9 which maps indications to categories


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

par(mar=c(7,7,3,1))
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
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


unecessary_message = dev.off()


n_active_total = sum(active_ti$cat %in% active_clinical$cat)
n_active_supported = sum(active_ti$cat %in% active_clinical$cat & active_ti$target_status=='genetically supported target')
n_hist_total = sum(hist_ti$cat %in% active_clinical$cat)
n_hist_supported = sum(hist_ti$cat %in% active_clinical$cat & hist_ti$target_status=='genetically supported target')

ctable = matrix(c(n_hist_total-n_hist_supported, n_hist_supported, n_active_total-n_active_supported, n_active_supported), nrow=2, byrow=T)
fisher_obj = fisher.test(ctable)

write(paste('T-I in Phase I-III with genetic support: ',n_hist_supported,'/',n_hist_total,' historical vs. ',
            n_active_supported,'/',n_active_total,' active. OR = ',
            formatC(fisher_obj$estimate,digits=2,format='fg'),', P = ',formatC(fisher_obj$p.value,digits=2,format='fg'),'\n',sep=''),text_stats_path,append=T)

####
# Figure 3
####

cat(file=stderr(), 'done.\nCreating Figure 3...')

resx=300
png(paste0(output_path,'/figure-3.png'),width=6.5*resx,height=9*resx,res=resx)

layout_matrix = matrix(c(1,1,
                         2,4,
                         3,4), nrow=3, byrow=T)
layout(layout_matrix, heights=c(2,.6,.6), widths=c(1,.9))
panel = 1

heavy_hitters = read_tsv('data/otg_heavy_hitter_curation.tsv', col_types=cols())
assoc$hh = assoc$original_trait %in% heavy_hitters$trait_reported

active_ti$gensup = active_ti$target_status == 'genetically supported target'
active_ti_forest = advancement_forest(active_ti,phase='active')
hist_ti_forest = advancement_forest(hist_ti,phase='historical')
historical_clinical_pg_nu = sum(hist_ti_forest$numerator[hist_ti_forest$label %in% active_clinical$cat])
historical_clinical_pg_de = sum(hist_ti_forest$denominator[hist_ti_forest$label %in% active_clinical$cat])
active_clinical_pg_nu = sum(active_ti_forest$numerator[active_ti_forest$label %in% active_clinical$cat])
active_clinical_pg_de = sum(active_ti_forest$denominator[active_ti_forest$label %in% active_clinical$cat])

pp$at_least_clinical = pp$ccat %in% c('Phase I','Phase II','Phase III','Launched')

indic %>% filter(genetic_insight != 'none') -> indic_insight
sim8 = sim[sim$comb_norm >= 0.8,]
assoc %>% filter(source != 'OTG' | l2g_share >= 0.5) -> assoc_l2g5
assoc_l2g5$associated_mesh_term = assoc_l2g5$mesh_term
sim8 %>%
  inner_join(assoc_l2g5, by=c('meshcode_a'='mesh_id')) %>%
  inner_join(indic_insight, by=c('meshcode_b'='indication_mesh_id')) %>%
  group_by(gene, assoc_mesh_id=meshcode_a, indication_mesh_id=meshcode_b, source) %>%
  summarize(.groups='keep', 
            assoc_mesh_term = min(associated_mesh_term),
            indication_mesh_term = min(indication_mesh_term),
            similarity = max(comb_norm),
            assoc_info = max(extra_info)) %>%
  ungroup() %>%
  mutate(ti_uid = paste0(gene,'-',indication_mesh_id)) -> all_possible_gensup_ti_assocs_all
all_possible_gensup_ti_assocs_all %>%
  group_by(ti_uid, gene, indication_mesh_id) %>%
  summarize(.groups='keep', 
            indication_mesh_term=min(indication_mesh_term), 
            sources=toString(unique(source))) %>%
  ungroup() %>%
  mutate(pp = ti_uid %in% pp$ti_uid[pp$at_least_clinical]) -> all_possible_gensup_ti_all


ctable = matrix(rep(as.integer(0),4), nrow=2, byrow=T, dimnames=list(gensup=c('0','1'), development=c('0','1')))

ctable['1','1'] = length(unique(pp$ti_uid[pp$at_least_clinical & pp$ti_uid %in% all_possible_gensup_ti_all$ti_uid]))
ctable['0','1'] = length(unique(pp$ti_uid[pp$at_least_clinical])) - ctable['1','1']
ctable['1','0'] = nrow(all_possible_gensup_ti_all) - ctable['1','1']
ctable['0','0'] = (nrow(universe)*sum(indic$genetic_insight!='none'))  - ctable['1','1'] -  - ctable['0','1'] -  - ctable['1','0']

gensup_dev_fisher_or = fisher.test(ctable)$estimate
gensup_dev_fisher_p = fisher.test(ctable)$p.value

write(paste("Fisher's test for enrichment of genetically supported T-I pairs among those developed: OR= ",
            formatC(gensup_dev_fisher_or,digits=1,format='f'),
            ', P = ',formatC(gensup_dev_fisher_p,digits=2,format='e'),'\n',sep=''),text_stats_path,append=T)

write(paste("Space of possible T-I pairs: ",nrow(universe)," genes times ",
            sum(indic$genetic_insight!='none')," indications with genetic insight = ",
            formatC((nrow(universe)*sum(indic$genetic_insight!='none')), big.mark=','),' T-I pairs','\n',sep=''),text_stats_path,append=T)

write(paste("Number of possible T-I pairs with genetic support: ",nrow(all_possible_gensup_ti_all),
            " = ",percent(nrow(all_possible_gensup_ti_all)/(nrow(universe)*sum(indic$genetic_insight!='none')),digits=2),'\n',sep=''),text_stats_path,append=T)

write(paste("Number of clinically developed T-I pairs: ",length(unique(pp$ti_uid[pp$at_least_clinical])),'\n',sep=''),text_stats_path,append=T)

write(paste("Number of clinically developed T-I pairs with genetic support: ",
            length(unique(pp$ti_uid[pp$at_least_clinical & pp$ti_uid %in% all_possible_gensup_ti_all$ti_uid])),
            " = ",percent(length(unique(pp$ti_uid[pp$at_least_clinical & pp$ti_uid %in% all_possible_gensup_ti_all$ti_uid]))/nrow(all_possible_gensup_ti_all),digits=2),
            " of possible / ",percent(length(unique(pp$ti_uid[pp$at_least_clinical & pp$ti_uid %in% all_possible_gensup_ti_all$ti_uid]))/length(unique(pp$ti_uid[pp$at_least_clinical])),digits=2),
            " of Phase I+",'\n',sep=''),text_stats_path,append=T)

# now restrict to *most* similar indication per each association
assoc %>% filter(source != 'OTG' | (l2g_share >= 0.5 & !hh))  -> assoc_l2g5_nohh
assoc_l2g5_nohh$associated_mesh_term = assoc_l2g5_nohh$mesh_term
sim8 %>%
  inner_join(assoc_l2g5_nohh, by=c('meshcode_a'='mesh_id')) %>%
  inner_join(indic_insight, by=c('meshcode_b'='indication_mesh_id')) %>%
  group_by(gene, assoc_mesh_id=meshcode_a, indication_mesh_id=meshcode_b, source) %>%
  summarize(.groups='keep', 
            assoc_mesh_term = min(associated_mesh_term),
            indication_mesh_term = min(indication_mesh_term),
            similarity = max(comb_norm),
            assoc_info = max(extra_info)) %>%
  ungroup() %>%
  arrange(gene, assoc_mesh_id, desc(similarity)) %>% # arrange by descending similarity within each gene-assoc
  group_by(gene, assoc_mesh_id) %>% 
  slice(1) %>% # take only the top most similar indication
  select(gene,
         assoc_mesh_id, 
         assoc_mesh_term, 
         indication_mesh_id, 
         indication_mesh_term, 
         similarity, 
         source, 
         assoc_info) %>%
  ungroup() %>%
  mutate(ti_uid = paste0(gene,'-',indication_mesh_id)) ->  all_possible_gensup_ti_assocs

all_possible_gensup_ti_assocs %>%
  group_by(ti_uid, gene, indication_mesh_id) %>%
  summarize(.groups='keep', term=min(indication_mesh_term), sources=toString(unique(source))) %>%
  ungroup() %>%
  mutate(pp = ti_uid %in% pp$ti_uid[pp$at_least_clinical]) -> all_possible_gensup_ti

write(paste("After removing heavy hitters and restricting to most similar indications, the number of possible T-I pairs with genetic support: ",nrow(all_possible_gensup_ti),
            " = ",percent(nrow(all_possible_gensup_ti)/(nrow(universe)*sum(indic$genetic_insight!='none')),digits=2),'\n',sep=''),text_stats_path,append=T)

write(paste("Number of clinically developed T-I pairs with genetic support: ",
            length(unique(pp$ti_uid[pp$at_least_clinical & pp$ti_uid %in% all_possible_gensup_ti$ti_uid])),
            " = ",percent(length(unique(pp$ti_uid[pp$at_least_clinical & pp$ti_uid %in% all_possible_gensup_ti$ti_uid]))/nrow(all_possible_gensup_ti),digits=2),'\n',sep=''),text_stats_path,append=T)


genelists = data.frame(x=c(-0.5,1:2,seq(3.5,7.5,1)),
                       list=c('all_genes','ab_tractable','sm_tractable','rhodop_gpcr','nuclear_receptors','enzymes','ion_channels','kinases'),
                       disp=c('all genes','predicted Ab tractable','predicted SM tractable','rhodopsin-like GPCRs','nuclear receptors','enzymes','ion channels','kinases'))

add_genelist_cols = function(tbl, genelistdf) {
  for (i in 1:nrow(genelistdf)) {
    if (substr(genelistdf$list[i],1,3)=='all') {
      tbl[,genelists$list[i]] = TRUE
    } else {
      genes = read.table(paste0('data/gene_lists/',genelists$list[i],'.tsv'),sep='\t',header=F)$V1
      tbl[,genelists$list[i]] = tbl$gene %in% genes
    }
  }
  return (tbl)
}

all_possible_gensup_ti = add_genelist_cols(all_possible_gensup_ti, genelists)

for (i in 1:nrow(areas)) {
  all_possible_gensup_ti[,areas$area[i]] = all_possible_gensup_ti$indication_mesh_id %in% indic_topl_match$indication_mesh_id[indic_topl_match$topl==areas$topl[i]]
}

all_possible_gensup_ti_assocs %>% filter(source=='intOGen') -> all_possible_gensup_ti_assocs_intogen

all_possible_gensup_ti_assocs_intogen %>%
  group_by(ti_uid, gene, indication_mesh_id) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() -> all_possible_gensup_ti_intogen

# note this table is already grouped by gene and the modal mechanism is chosen:
intogen_mechanisms = read_tsv('data/intogen_genes.tsv', col_types=cols())

# otherwise you get oddities like TP53 is a tumor suppressor in most histologic types but due to some
# random error was classified as oncogene in one tissue so it turns up in the oncogene list
all_possible_gensup_ti_intogen$direction = intogen_mechanisms$mechanism[match(all_possible_gensup_ti_intogen$gene, intogen_mechanisms$gene)]
all_possible_gensup_ti_intogen$direction[is.na(all_possible_gensup_ti_intogen$direction)] = 'unknown'
all_possible_gensup_ti_intogen$pp = all_possible_gensup_ti_intogen$ti_uid %in% pp$ti_uid[pp$at_least_clinical]

all_possible_gensup_ti_intogen = add_genelist_cols(all_possible_gensup_ti_intogen, genelists)

intogen_directions = data.frame(dir=c('oncogene','tumor_suppressor','unknown'),y=3:1,disp=c('oncogene','tumor suppressor','unknown'))

for (i in 1:nrow(intogen_directions)) {
  all_possible_gensup_ti_intogen[,intogen_directions$dir[i]] = all_possible_gensup_ti_intogen$direction == intogen_directions$disp[i]
}

for (i in 1:nrow(areas)) {
  all_possible_gensup_ti_intogen[,areas$area[i]] = all_possible_gensup_ti_intogen$indication_mesh_id %in% indic_topl_match$indication_mesh_id[indic_topl_match$topl==areas$topl[i]]
}

ti_utilization = function(posstable, geneset, indicset) {
  denominator = sum(posstable[,geneset] & posstable[,indicset])
  numerator = sum(posstable[,geneset] & posstable[,indicset] & posstable$pp)
  return (c(numerator, denominator))
}

target_utilization = function(posstable, geneset, indicset) {
  denominator = length(unique(posstable$gene[posstable[,geneset] & posstable[,indicset]]))
  numerator = length(unique(posstable$gene[posstable[,geneset] & posstable[,indicset] & posstable$pp]))
  return (c(numerator, denominator))
}

all_possible_gensup_ti$all_areas = TRUE
all_possible_gensup_ti_intogen$all_areas = TRUE

areas %>%
  select(area, y) %>%
  add_row(area='all_areas',y=nrow(areas)+1.5, .before=1) %>%
  mutate(disp = gsub('_',' ',area)) -> area_meta

panel = 1

maxcol = '#6E016B'
ylims = range(area_meta$y) + c(-0.5, 0.5)
xlims = range(genelists$x) + c(-0.5, 0.5)
par(mar=c(0.5,8,8,2))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=2, at=area_meta$y, labels=area_meta$disp, lwd=0, line=-0.25, las=2, cex.axis=1.1)
par(xpd=T)
text(x=genelists$x, y=rep(max(ylims), nrow(genelists))+0.1, labels=genelists$disp, adj=c(0,0), srt=45, cex=1.1)
par(xpd=F)
utilization_table = tibble(area=character(0), 
                           genelist=character(0), 
                           developed=integer(0),
                           supported=integer(0),
                           utilization=numeric(0))
for (a in 1:nrow(area_meta)) {
  for (g in 1:nrow(genelists)) {
    numden = ti_utilization(all_possible_gensup_ti, genelists$list[g], area_meta$area[a])
    utilization = numden[1] / numden[2]
    if (numden[2] == 0) utilization = 0
    utilization_table = rbind(utilization_table,
                              tibble(area=area_meta$area[a],
                                     genelist=genelists$disp[g],
                                     developed=numden[1],
                                     supported=numden[2],
                                     utilization=utilization))
    disp = gsub(' ','',paste0(formatC(utilization*100, format='fg', digits=2),'%'))
    col = alpha(maxcol, utilization)
    rect(xleft=genelists$x[g]-0.5, xright=genelists$x[g]+0.5, ybottom=area_meta$y[a]-0.5, ytop=area_meta$y[a]+0.5, col=col, border=NA)
    text(x=genelists$x[g], y=area_meta$y[a]-0.2, labels=disp, cex=1.2, font=2)
    text(x=genelists$x[g], y=area_meta$y[a]-0.1, labels=paste0(numden[1],'/',formatC(numden[2],format='fg',big.mark=',')), pos=3, cex=0.8)
  }
}
abline(h=mean(area_meta$y[1:2]),lwd=0.5)
abline(v=0.25,lwd=0.5)
abline(v=2.75,lwd=0.5)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(utilization_table, 'Proportion of all possible genetically supported target-indication pairs that have been developed, by therapy area and gene list.')

maxcol = '#6E016B'
ylims = range(intogen_directions$y) + c(-0.5, 0.5)
xlims = range(genelists$x) + c(-0.5, 0.5)
par(mar=c(1,11,7,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=2, at=intogen_directions$y, labels=intogen_directions$disp, lwd=0, las=2, line=-0.25, cex.axis=1.1)
par(xpd=T)
text(x=genelists$x, y=rep(max(ylims), nrow(genelists))+0.05, labels=genelists$disp, adj=c(0,0), srt=45, cex=0.9)
par(xpd=F)
intogen_utilization_table = tibble(area=character(0), 
                           genelist=character(0), 
                           developed=integer(0),
                           supported=integer(0),
                           utilization=numeric(0))
for (a in 1:nrow(intogen_directions)) {
  for (g in 1:nrow(genelists)) {
    numden = target_utilization(all_possible_gensup_ti_intogen, genelists$list[g], intogen_directions$dir[a])
    utilization = numden[1] / numden[2]
    if (numden[2]==0) { # denominator 0, undefined
      disp = '-'
      col = alpha(maxcol, 0)
      utilization = 0
    } else {
      disp = gsub(' ','',paste0(formatC(utilization*100, format='fg', digits=2),'%'))
      col = alpha(maxcol, utilization)
    }
    intogen_utilization_table = rbind(intogen_utilization_table,
                              tibble(area=intogen_directions$dir[a],
                                     genelist=genelists$disp[g],
                                     developed=numden[1],
                                     supported=numden[2],
                                     utilization=utilization))
    rect(xleft=genelists$x[g]-0.5, xright=genelists$x[g]+0.5, ybottom=intogen_directions$y[a]-0.5, ytop=intogen_directions$y[a]+0.5, col=col, border=NA)
    text(x=genelists$x[g], y=intogen_directions$y[a]-0.2, labels=disp, cex=0.8, font=2)
    text(x=genelists$x[g], y=intogen_directions$y[a]-0.15, labels=paste0(numden[1],'/',numden[2]), pos=3, cex=0.6)
  }
}
mtext(side=2, line=8.5, text='IntOGen\nmechanism', cex=0.7)
abline(v=0.25,lwd=0.5)
abline(v=2.75,lwd=0.5)
mtext(letters[panel], side=3, cex=2, adj = -0.3, line = 0.5)
panel = panel + 1

write_supp_table(intogen_utilization_table, 'Proportion of all possible IntOGen-supported target-indication pairs that have been developed, by driver mechanism and gene list.')



## gsup_bins setup
pp$ccatnum = meta_ccat$num[match(pp$ccat, meta_ccat$cat)]
pp$ti_uid = paste0(pp$gene,'-',pp$indication_mesh_id)
all_possible_gensup_ti$ti_uid = paste0(all_possible_gensup_ti$gene,'-',all_possible_gensup_ti$indication_mesh_id)
all_possible_gensup_ti$gene_is_pp_clinical = all_possible_gensup_ti$gene %in% pp$gene[pp$at_least_clinical]
all_possible_gensup_ti %>%
  filter(gene_is_pp_clinical) %>%
  group_by(gene) %>%
  summarize(.groups='keep',
            n_gsup_ti = n(), 
            n_gsup_pp_ti = sum(pp)) -> gsup_targets
pp$gene_is_gsup_target = pp$gene %in% all_possible_gensup_ti$gene
pp$ti_is_nonsup = !(pp$ti_uid %in% all_possible_gensup_ti$ti_uid)
pp %>%
  filter(gene_is_gsup_target & ti_is_nonsup & at_least_clinical) %>%
  group_by(gene) %>%
  summarize(.groups='keep', n_nonsup_pp_ti = length(unique(ti_uid))) -> nonsup_targets_pp
gsup_targets$n_nonsup_pp_ti = nonsup_targets_pp$n_nonsup_pp_ti[match(gsup_targets$gene, nonsup_targets_pp$gene)]
gsup_targets$n_nonsup_pp_ti[is.na(gsup_targets$n_nonsup_pp_ti)] = 0
gsup_targets$launched = gsup_targets$gene %in% pp$gene[pp$ccat=='Launched']
gsup_targets$omim = gsup_targets$gene %in% assoc$gene[assoc$source=='OMIM']
gsup_targets$gwas = gsup_targets$gene %in% assoc$gene[assoc$source=='OTG' & assoc$l2g_share >= 0.5 & !is.na(assoc$l2g_share)]

write_supp_table(gsup_targets, 'Number of possible genetically supported indications, supported developed indications, and unsupported developed indications, by target.')

write(paste('Of targets with both ≥1 supported indication and ≥1 Phase I+ indication, ',
            percent(mean(gsup_targets$n_gsup_pp_ti==0)),' (',sum(gsup_targets$n_gsup_pp_ti==0),'/',nrow(gsup_targets),
            ') have been developed only for unsupported indications.'
            ,'\n',sep=''),text_stats_path,append=T)

write(paste('Of targets with both ≥1 supported indication and ≥1 Phase I+ indication, ',
            percent(mean(gsup_targets$n_nonsup_pp_ti==0)),' (',sum(gsup_targets$n_nonsup_pp_ti==0),'/',nrow(gsup_targets),
            ') have been developed only for supported indications.'
            ,'\n',sep=''),text_stats_path,append=T)

write(paste('Of targets with both ≥1 supported indication and ≥1 Phase I+ indication, ',
            percent(mean(gsup_targets$n_nonsup_pp_ti > 0 & gsup_targets$n_gsup_pp_ti > 0)),' (',sum(gsup_targets$n_nonsup_pp_ti > 0 & gsup_targets$n_gsup_pp_ti > 0),'/',nrow(gsup_targets),
            ') have been developed for both.'
            ,'\n',sep=''),text_stats_path,append=T)

write(paste('Of targets with both ≥1 supported indication and ≥1 Phase I+ indication, there were ',
            sum(gsup_targets$n_nonsup_pp_ti),' unsupported indications and ',sum(gsup_targets$n_gsup_pp_ti),' supported indications pursued, a ratio of ',
            formatC(sum(gsup_targets$n_nonsup_pp_ti)/sum(gsup_targets$n_gsup_pp_ti), format='f', digits=1),'\n',sep=''),text_stats_path,append=T)

write(paste('Of launched targets with both ≥1 supported indication and ≥1 Phase I+ indication, there were ',
            sum(gsup_targets$n_nonsup_pp_ti[gsup_targets$launched]),' unsupported indications and ',sum(gsup_targets$n_gsup_pp_ti[gsup_targets$launched]),' supported indications pursued, a ratio of ',
            formatC(sum(gsup_targets$n_nonsup_pp_ti[gsup_targets$launched])/sum(gsup_targets$n_gsup_pp_ti[gsup_targets$launched]), format='f', digits=1),'\n',sep=''),text_stats_path,append=T)

write(paste('Of non-launched targets with both ≥1 supported indication and ≥1 Phase I+ indication, there were ',
            sum(gsup_targets$n_nonsup_pp_ti[!gsup_targets$launched]),' unsupported indications and ',sum(gsup_targets$n_gsup_pp_ti[!gsup_targets$launched]),' supported indications pursued, a ratio of ',
            formatC(sum(gsup_targets$n_nonsup_pp_ti[!gsup_targets$launched])/sum(gsup_targets$n_gsup_pp_ti[!gsup_targets$launched]), format='f', digits=1),'\n',sep=''),text_stats_path,append=T)

ctable = matrix(c(sum(gsup_targets$n_nonsup_pp_ti[!gsup_targets$launched]),
                  sum(gsup_targets$n_gsup_pp_ti[!gsup_targets$launched]),
                  sum(gsup_targets$n_nonsup_pp_ti[gsup_targets$launched]),
                  sum(gsup_targets$n_gsup_pp_ti[gsup_targets$launched])), nrow=2, byrow=T)
fisher_obj = fisher.test(ctable)

write(paste('Fisher test for difference between these ratios: OR=',
            formatC(fisher_obj$estimate, format='f', digits=1),', P=',formatC(fisher_obj$p.value, format='f', digits=2),'.',
            '\n',sep=''),text_stats_path,append=T)



# View(gsup_targets[gsup_targets$n_gsup_pp_ti==0 & !gsup_targets$launched & gsup_targets$gwas & gsup_targets$omim,])
all_possible_gensup_ti$ti_is_pp_clinical = all_possible_gensup_ti$ti_uid %in% pp$ti_uid[pp$at_least_clinical]
all_possible_gensup_ti %>%
  filter(gene_is_pp_clinical) %>%
  group_by(ti_uid, gene, indication_mesh_id) %>%
  summarize(.groups='keep', pp=ifelse(sum(ti_is_pp_clinical) > 0, 1, 0)) -> gsup_pp_ti

maxbin = 10
gsup_targets$bin_n_gsup_ti = pmin(gsup_targets$n_gsup_ti,maxbin)
gsup_targets %>%
  group_by(bin_n_gsup_ti) %>%
  summarize(.groups='keep',
            mean_sup = mean(n_gsup_pp_ti),
            sd_sup = sd(n_gsup_pp_ti),
            mean_non = mean(n_nonsup_pp_ti),
            sd_non = sd(n_nonsup_pp_ti),
            n=n()) -> gsup_bins
gsup_bins$disp = as.character(gsup_bins$bin_n_gsup_ti)
gsup_bins$disp[gsup_bins$bin_n_gsup_ti==maxbin] = paste0(as.character(maxbin),'+')
gsup_bins$l95_sup = gsup_bins$mean_sup - 1.96 * gsup_bins$sd_sup/sqrt(gsup_bins$n)
gsup_bins$u95_sup = gsup_bins$mean_sup + 1.96 * gsup_bins$sd_sup/sqrt(gsup_bins$n)
gsup_bins$l95_non = gsup_bins$mean_non - 1.96 * gsup_bins$sd_non/sqrt(gsup_bins$n)
gsup_bins$u95_non = gsup_bins$mean_non + 1.96 * gsup_bins$sd_non/sqrt(gsup_bins$n)

supcol = '#43A2CA'
noncol = '#BC7642'
binwidth = 0.4
xlims=c(-20,20)
ylims=c(0,maxbin+0.5)
par(mar=c(3.5,4,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=-20:20, labels=NA, tck=-0.025)
axis(side=1, at=c(-20, -10, 10, 20), labels=NA, tck=-0.05)
axis(side=1, at=c(-20, -10, 10, 20), labels=c(20, 10, 10, 20), lwd=0, line=-0.75)
axis(side=1, at=0, labels=NA, tck=-0.15)
mtext(side=1, line=1.25, at=c(-10, 10), text=c('Supported','Unsupported'), cex=0.8, col=c(supcol, noncol))
mtext(side=1, line=2.25, at=0, text=c('Indications pursued to phase I+'), cex=0.8)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
axis(side=2, at=gsup_bins$bin_n_gsup_ti, tck=-0.02, labels=NA)
axis(side=2, at=c(0,5,10), labels=NA, tck=-0.05)
axis(side=2, at=c(0,5,10), labels=c(0,5,10), line=-0.5, lwd=0, las=2)
mtext(side=2, line=1.5, text='Possible indications\nwith genetic support', cex=0.8)
abline(v=0)
rect(xleft=-gsup_bins$mean_sup, xright=0, ybottom=gsup_bins$bin_n_gsup_ti-binwidth, ytop=gsup_bins$bin_n_gsup_ti+binwidth, col=supcol, border=NA)
arrows(x0=-gsup_bins$u95_sup, x1=-gsup_bins$l95_sup, y0=gsup_bins$bin_n_gsup_ti, y1=gsup_bins$bin_n_gsup_ti, angle=90, code=3, length=0.03)
rect(xleft=0, xright=gsup_bins$mean_non, ybottom=gsup_bins$bin_n_gsup_ti-binwidth, ytop=gsup_bins$bin_n_gsup_ti+binwidth, col=noncol, border=NA)
arrows(x0=gsup_bins$l95_non, x1=gsup_bins$u95_non, y0=gsup_bins$bin_n_gsup_ti, y1=gsup_bins$bin_n_gsup_ti, angle=90, code=3, length=0.03)
mtext(letters[panel], side=3, cex=2, adj = 0.0, line = 0.5)
panel = panel + 1

write_supp_table(gsup_bins, 'Histogram of supported vs. unsupported indications pursued for developed targets.')

hist_ti$gensup = hist_ti$target_status == 'genetically supported target'
drug_phase_summary$gensup = hist_ti$gensup[match(drug_phase_summary$ti_uid, hist_ti$ti_uid)]
drug_phase_summary = drug_phase_summary[!is.na(drug_phase_summary$gensup),]
drug_phase_summary %>%
  filter(maxphase != 'Preclinical') %>%
  group_by(phasenum, phase, maxphasenum, maxphase) %>%
  summarize(.groups='keep',
            n_gensup = sum(n_drugs[gensup]),
            n_total = sum(n_drugs)) %>%
  ungroup() %>%
  arrange(maxphasenum, phasenum) %>%
  mutate(y = max(row_number()) - row_number() + 1) -> d_by_max

d_by_max[,c('pg','pg_l95','pg_u95')] = binom.confint(d_by_max$n_gensup, d_by_max$n_total, method='wilson')[,c('mean','lower','upper')]

d_by_max %>%
  group_by(maxphase) %>%
  summarize(.groups='keep', miny=min(y), maxy=max(y), midy=mean(y)) %>%
  ungroup() %>%
  mutate(disp = gsub('clinical','',gsub('Phase ','',maxphase))) -> tranches
tranche_breaks = tranches$maxy + 0.5

xlims = c(0, 0.16)
ylims = range(d_by_max$y) + c(-0.5, 0.51)
par(mar=c(3,9,4,6))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=1, at=0:100/100, labels=NA, tck=-0.025)
axis(side=1, at=0:20/20, labels=NA, tck=-0.05)
axis(side=1, at=0:20/20, labels=percent(0:20/20,digits=0), lwd=0, line=-0.5)
mtext(side=1, line=1.6, cex=1, text='P(G)')
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
outer_line = 5
overhang = 0.4
for (i in 1:nrow(tranches)) {
  axis(side=2, line=outer_line, tck=0.02, at=c(tranches$miny[i]-overhang, tranches$maxy[i]+overhang), labels=NA)
}
mtext(side=2, line=outer_line+0.25, at=tranches$midy, text=tranches$maxphase, cex=.6, las=2)
mtext(side=2, line=0.5, at=d_by_max$y, text=d_by_max$phase, cex=.6, las=2)
mtext(side=2, at=max(ylims)+1.5, line=outer_line + 0.5, las=2, text='Drug\nhighest\nphase\nreached',font=2,cex=0.6)
mtext(side=2, at=max(ylims)+1.5, line=0.5, las=2, text='Indication\nphase', font=2, cex=0.6)
mtext(side=4, at=max(ylims)+1.5, line=0.5, las=2, text='Supported/\ntotal', font=2, cex=0.6)
mtext(side=4, at=d_by_max$y, text=paste0(formatC(d_by_max$n_gensup,big.mark=','),'/',formatC(d_by_max$n_gensup+d_by_max$n_total,format='d',big.mark=',')), cex=0.6, las=2, line=0.25)
abline(h=tranche_breaks, lwd=0.25)
points(d_by_max$pg, d_by_max$y, pch=19) # means
segments(x0=d_by_max$pg_l95, x1=d_by_max$pg_u95, y0=d_by_max$y, lwd=2) # 95%CIs

mtext(letters[panel], side=3, cex=2, adj = 0.0, line = 0.5)
panel = panel + 1

write_supp_table(d_by_max, "Proportion of indications with genetic support, by drug's highest phase reached and phase for each indication")

as_tibble(lapply(drug_phase_summary, rep, drug_phase_summary$n_drugs)) %>%
  select(phase, maxphase, gensup) %>%
  mutate(phase = factor(phase, ordered=T, levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'))) %>%
  mutate(maxphase = factor(maxphase, ordered=T, levels=c('Preclinical','Phase I','Phase II','Phase III','Launched'))) -> ordinal_data

ordinal_model = polr(phase ~ maxphase + gensup, data=ordinal_data, Hess=T)
ordinal_model_coefs = coef(summary(ordinal_model))
ordinal_model_betas = ordinal_model_coefs[,'Value']
ordinal_model_p = pnorm(abs(ordinal_model_coefs[,'t value']), lower.tail = FALSE) * 2
ordinal_repur_model_smry = tibble(variable=rownames(ordinal_model_coefs),
                                  beta=ordinal_model_betas, 
                                  p=ordinal_model_p)

write(paste('Ordinal logistic model for drug-indication advancement: genetic support beta= ',
            formatC(ordinal_repur_model_smry$beta[ordinal_repur_model_smry$variable=='gensupTRUE'],digits=2,format='fg'),
            ', P = ',
            formatC(ordinal_repur_model_smry$p[ordinal_repur_model_smry$variable=='gensupTRUE'],digits=2,format='e'),'\n',sep=''),
      text_stats_path,append=T)


unecessary_message = dev.off()



####
# FIGURE S7
####

cat(file=stderr(), 'done.\nCreating Figure S7...')

resx=300
png(paste0(output_path,'/figure-s7.png'),width=6.5*resx,height=8*resx,res=resx)

layout_matrix = matrix(c(1,2,3,3),nrow=2,byrow=T)
layout(layout_matrix, heights=c(1,3))

panel = 1

gtcolor = '#3489CA'
ticolor = '#CB7635'

par(mar=c(3,3,3,1))
xlims = c(2007, 2022)
xats = 2007:2022
xbigs = c(2010, 2015, 2020)
ylims = c(0, 81000)
yats = 0:10*1e4
ybigs = c(0, 20000, 40000, 60000, 80000)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, line=-0.5, lwd=0, labels=xbigs)
mtext(side=1, line=1.5, text='Year')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, line=-0.5, lwd=0, labels=paste0(ybigs/1000,'K'), las=2)
mtext(side=2, line=2, text='Cumulative pairs')
points(gtyc$min_year, gtyc$cumn, col=gtcolor, lwd=1.5, type='l')
points(tiyc$min_year, tiyc$cumn, col=ticolor, lwd=1.5, type='l')
legend('topleft', c('G-I pairs supported','G-T associations'), col=c(ticolor,gtcolor), lwd=1.5, bty='n', cex=1)
mtext(letters[panel], side=3, cex=2, adj = 0.0, line = 0.5)
panel = panel + 1

plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')


maxcol = '#6E016B'
ylims = range(area_meta$y) + c(-0.5, 0.5)
xlims = range(genelists$x) + c(-0.5, 0.5)
par(mar=c(1,10,6,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=2, at=area_meta$y, labels=area_meta$disp, lwd=0, las=2)
par(xpd=T)
text(x=genelists$x, y=rep(max(ylims), nrow(genelists)), labels=genelists$disp, adj=c(0,0), srt=45)
par(xpd=F)
target_utilization_table = tibble(area=character(0), 
                                  genelist=character(0), 
                                  developed=integer(0),
                                  supported=integer(0),
                                  utilization=numeric(0))
for (a in 1:nrow(area_meta)) {
  for (g in 1:nrow(genelists)) {
    numden = target_utilization(all_possible_gensup_ti, genelists$list[g], area_meta$area[a])
    utilization = numden[1] / numden[2]
    if (numden[2]==0) utilization = 0
    target_utilization_table = rbind(target_utilization_table,
                                     tibble(area=area_meta$area[a],
                                            genelist=genelists$disp[g],
                                            developed=numden[1],
                                            supported=numden[2],
                                            utilization=utilization))
    disp = gsub(' ','',paste0(formatC(utilization*100, format='fg', digits=2),'%'))
    col = alpha(maxcol, utilization)
    rect(xleft=genelists$x[g]-0.5, xright=genelists$x[g]+0.5, ybottom=area_meta$y[a]-0.5, ytop=area_meta$y[a]+0.5, col=col, border=NA)
    text(x=genelists$x[g], y=area_meta$y[a]-0.2, labels=disp, cex=1, font=2)
    text(x=genelists$x[g], y=area_meta$y[a]-0.1, labels=paste0(numden[1],'/',numden[2]), pos=3, cex=0.7)
  }
}
abline(h=mean(area_meta$y[1:2]),lwd=0.5)
abline(v=c(0.25,2.75),lwd=0.5)
mtext(letters[panel], side=3, cex=2, adj = 0.0, line = 0.5)

write_supp_table(target_utilization_table, 'Proportion of all possible genetically supported targets that have been developed, by therapy area and gene list.')

unecessary_message = dev.off()

cat(file=stderr(), 'done.\nFinalizing supplementary tables...')

#####
# SUPPLEMENT
#####

# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)



elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))
