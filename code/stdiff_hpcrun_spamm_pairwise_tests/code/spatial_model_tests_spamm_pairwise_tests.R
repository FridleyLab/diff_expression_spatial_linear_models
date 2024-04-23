# Helpers ----------------------------------------------------------------------

##
# non_spatial_de_hpc
#
non_spatial_de_hpc = function(expr_data=NULL, annot1=NULL, annot2=NULL){
  # Recode clusters
  meta1_tmp = annot1
  meta2_tmp = annot2
  expr_tmp = expr_data %>%
    dplyr::select(group, meta, xpos, ypos, exprval) %>%
    dplyr::filter(meta %in% c(meta1_tmp, meta2_tmp))
    # dplyr::mutate(meta=factor(dplyr::case_when(meta != meta_tmp ~ 'other',
    #                                            TRUE ~ meta), levels=c('other', meta_tmp)))
  
  # Create non-spatial model
  start_nonsp_t = Sys.time()
  res_ls = spaMM::fitme(formula=as.formula('exprval~meta'), data=expr_tmp, method="REML")
  end_nonsp_t = difftime(Sys.time(), start_nonsp_t, units='min')
  
  return(list(nonspmod=res_ls,
              exectime=end_nonsp_t))
}


##
# spatial_de_hpc
#
spatial_de_hpc = function(non_sp_mods=NULL){
  # Initiates list to store results
  res_ls = list()
  
  # Extract data from the same model object
  expr_tmp = non_sp_mods[['data']]
  
  # Update model call to use nlme package in new environemnt
  #non_sp_mods[['call']] = str2lang('nlme::lme.formula(fixed=as.formula(paste0("exprval~meta")), data=expr_tmp, random=~1|group, method="REML")')
  
  # Run spherical model and catch error if no convergence
  # start_sph_t = Sys.time()
  # sph_out = tryCatch({
  #   update(non_sp_mods, correlation=nlme::corSpher(form=~xpos+ypos, nugget=T, metric='euclidean'), method='REML', control=nlme::lmeControl(returnObject=F, msVerbose=F, gradHess=T))
  # }, error=function(err){return(err)})
  # end_sph_t = difftime(Sys.time(), start_sph_t, units='min')
  # 
  # if(any(class(sph_out) == 'simpleError')){
  #   res_ls[['sph']][[1]] = 'no_conv'
  #   res_ls[['sph']][[2]] = names(non_sp_mods$coefficients$fixed)
  # } else{
  #   res_ls[['sph']] = sph_out
  # }
  
  # Run exponential model and catch error if no convergence
  start_exp_t = Sys.time()
  # Run exponential model and catch error if no convergence
  exp_out = tryCatch({
    spaMM::fitme(formula=as.formula(paste0("exprval~meta+Matern(1|xpos+ypos)")), 
                 data=expr_tmp, fixed=list(nu=0.5), method="REML", 
                 control.HLfit=list(algebra="decorr"))  }, error=function(err){return(err)})
  end_exp_t = difftime(Sys.time(), start_exp_t, units='min')
  
  if(any(class(exp_out) == 'simpleError')){
    res_ls[['exp']][[1]] = 'no_conv'
    res_ls[['exp']][[2]] = names(non_sp_mods$coefficients$fixed)
  } else{
    res_ls[['exp']] = exp_out
  }
  
  return(list(spmod=res_ls,
              exectime=list(#sph_time=end_sph_t,
                            exp_time=end_exp_t)))
}

   
##
# ANALYSIS BEGINS 
# Differential gene expression of clusters using non-spatial models and spatial 
# models with spatial covariance structures
#
library('magrittr')
library('spatialGE')
#devtools::load_all('../../../spatialGE/')

cmdargs = commandArgs(trailingOnly=T)
genename = as.character(cmdargs[1])
#genename = 'ABCA1'
cluster1 = as.character(cmdargs[2])
cluster2 = as.character(cmdargs[3])
#cluster = 'none'
stlist_fp = as.character(cmdargs[4])
#stlist_fp = '../data/visium_stlist_w_clusters_raviglioblastoma.RDS'
samplename = as.character(cmdargs[5])
#samplename = 'UKF243'

# Load STlist
x = readRDS(stlist_fp)

# Test sample is in STlist
if(length(grep(samplename, names(x@spatial_meta))) != 1){
  stop('The requested samples are not present in the STlist, or requested more than one sample')
}

# Test that annotation category is present in STlist
# if(!(any(x@spatial_meta[[samplename]][['annots']] == cluster))){
#   stop('The requested annotations are not present in the STlist, or requested more than one annotation')
# }

# Test that gene is present in STlist
if(!(any(rownames(x@tr_counts[[samplename]]) == genename))){
  stop('The requested annotations are not present in the STlist, or requested more than one annotation')
}

cat(crayon::yellow(paste0('Sample: ', samplename, '...\n')))
cat(crayon::yellow(paste0('Gene: ', genename, '...\n')))

# Create data frame with expression and coordinate data
# Add group dummy column and select relevant columns
# Add character to annotation column that it's treated a factor in case it's an integer
expr_df = as.data.frame(as.matrix(x@tr_counts[[samplename]])) %>%
  tibble::rownames_to_column(var='gene') %>%
  dplyr::filter(gene %in% genename) %>%
  tibble::column_to_rownames(var='gene') %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var='libname') %>%
  dplyr::right_join(x@spatial_meta[[samplename]] %>%
                      tibble::add_column(group=1, .after='libname') %>%
                      dplyr::select(libname, group, ypos, xpos, meta_orig=annots),. , by='libname') %>%
  tibble::column_to_rownames(var='libname') %>%
  dplyr::rename(exprval := !!genename)

# Create "dictionary" with coded annotations (to avoid potentially problematic characters)
meta_dict = tibble::tibble(orig_annot=unique(expr_df[['meta_orig']]),
                           coded_annot=paste0('c', 1:length(unique(expr_df[['meta_orig']]))))

# Recode annotations in expression data frame
for(spotrow in 1:nrow(expr_df)){
  expr_df[['meta']][spotrow] = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == expr_df[['meta_orig']][spotrow] ]
}
rm(spotrow) # Clean env
expr_df = expr_df %>% dplyr::select(-c('meta_orig')) # Remove original annotation and keep 'meta' as the coded annotations column

cat('Running non-spatial test...\n')
# Run models in parallel and get DE results
non_sp_models = non_spatial_de_hpc(expr_data=expr_df, 
                                   annot1=meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == cluster1 ],
                                   annot2=meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == cluster2 ]) # NEED TO MAKE FUNCTION TO SAVE EXEC TIME

# Calculate log-fold changes
meta1_tmp = unique(non_sp_models[['nonspmod']][['data']][['meta']]) %>% grep(meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == cluster2 ], ., value=T, invert=T)
meta2_tmp = unique(non_sp_models[['nonspmod']][['data']][['meta']]) %>% grep(meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == cluster1 ], ., value=T, invert=T)
avgexpr_ref = mean(non_sp_models[['nonspmod']][['data']] %>% dplyr::filter(meta != meta2_tmp) %>% dplyr::select(exprval) %>% unlist())
avgexpr_other = mean(non_sp_models[['nonspmod']][['data']] %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select(exprval) %>% unlist())

# Add results and logFC to list of data frames
non_sp_de = tibble::tibble(gene=genename,
                           cluster1=meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta1_tmp ],
                           cluster2=meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta2_tmp ],
                           avc_log2fc=(avgexpr_ref - avgexpr_other),
                           p_val=spaMM::summary.HLfit(non_sp_models[['nonspmod']], details=c(p_value="Wald"), verbose=F)[['beta_table']][2, 4],
                           nonsp_time_min=as.vector(non_sp_models[['exectime']])) %>% 
  tibble::add_column(sample=samplename, .before=1)

rm(meta1_tmp, meta2_tmp, avgexpr_ref, avgexpr_other) # Clean environment

cat('Running spatial test...\n')
# Run models
sp_models = spatial_de_hpc(non_sp_mods=non_sp_models[['nonspmod']])

rm(expr_df) # Clean environment

# Compile results of spatial models and merge with non-spatial results
# Check that model results exists and extract values, gene names, and clusters tested
# if(class(sp_models[['spmod']][['sph']]) == 'lme'){
#   sph_summ = summary(sp_models[['spmod']][['sph']])
#   sph_lrt = anova(non_sp_models[['nonspmod']], sp_models[['spmod']][['sph']])[['p-value']][2]
#   meta_sph_tmp = unique(sp_models[['spmod']][['sph']][['data']][['meta']]) %>% grep('other', ., value=T, invert=T)
#   meta_sph_tmp = meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta_sph_tmp ]
#   sph_time_tmp=as.vector(sp_models[['exectime']][['sph_time']])
# } else{
#   sph_summ = list(AIC=-9999, BIC=-9999, tTable=data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999)))
#   meta_sph_tmp = cluster
#   sph_time_tmp=NA_real_
#   sph_lrt=NA_real_
# }

# exp_lrt_obj = tryCatch({
#   lrt_in = spaMM::get_RLRsim_args(sp_models[['spmod']][['exp']], non_sp_models[['nonspmod']], verbose=F)
#   sim_lrt = RLRsim::LRTSim(lrt_in[['X']], lrt_in[['Z']], 1, lrt_in[["sqrt.Sigma"]], nsim=100)
#   obs_lrt = 2*(logLik(sp_models[['spmod']][['exp']])-logLik(non_sp_models[['nonspmod']]))
#   lrt_res = list(obs_lrt=as.vector(obs_lrt),
#                  lrt_pval=(sum(sim_lrt >= obs_lrt) + 1) / (length(sim_lrt) + 1))
# }, error=function(err){return(err)})

exp_lrt_obj = tryCatch({
  list(aic_nonsp=spaMM::extractAIC.HLfit(non_sp_models[['nonspmod']])[[2]],
       aic_sp=spaMM::extractAIC.HLfit(sp_models[['spmod']][['exp']])[[2]])
}, error=function(err){return(err)})

if(any(class(exp_lrt_obj) == 'simpleError')){
  exp_lrt_obj = list(aic_nonsp=NA_real_, aic_sp=NA_real_)
} #else{
#  exp_lrt_obj = exp_lrt_obj[['lrt_pval']]
#}

if(any(class(sp_models[['spmod']][['exp']]) == 'HLfit')){
  exp_summ = spaMM::summary.HLfit(sp_models[['spmod']][['exp']], details=c(p_value='Wald'), verbose=F)
  #exp_lrt = exp_lrt_obj
  non_sp_aic = exp_lrt_obj[['aic_nonsp']]
  sp_aic = exp_lrt_obj[['aic_sp']]
  meta1_exp_tmp = unique(sp_models[['spmod']][['exp']][['data']][['meta']]) %>% grep(meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == cluster2 ], ., value=T, invert=T)
  meta1_exp_tmp = meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta1_exp_tmp ]
  meta2_exp_tmp = unique(sp_models[['spmod']][['exp']][['data']][['meta']]) %>% grep(meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == cluster1 ], ., value=T, invert=T)
  meta2_exp_tmp = meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta2_exp_tmp ]
  exp_time_tmp=as.vector(sp_models[['exectime']][['exp_time']])
} else{
  exp_summ = list(#AIC=-9999, #BIC=-9999, 
                  beta_table=data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999)))
  #exp_lrt=NA_real_
  meta1_exp_tmp = cluster1
  meta2_exp_tmp = cluster2
  exp_time_tmp=NA_real_
  non_sp_aic = NA_real_
  sp_aic = NA_real_
}

# Make sure genes and clusters from both exponential and spherical models match
#if((meta_exp_tmp == meta_sph_tmp)){
  sp_de = tibble::tibble(sample=samplename,
                         gene=genename,
                         cluster1=meta1_exp_tmp,
                         cluster2=meta2_exp_tmp,
                         #sph_p_val=sph_summ[['tTable']][2, 5],
                         #sph_aic=sph_summ[['AIC']],
                         #sph_bic=sph_summ[['BIC']],
                         #sph_lrt_p_val=sph_lrt,
                         #sph_time_min=sph_time_tmp,
                         exp_p_val=exp_summ[['beta_table']][2, 4],
                         non_sp_aic=non_sp_aic,
                         sp_aic=sp_aic,
                         #exp_bic=exp_summ[['BIC']],
                         #exp_lrt_p_val=exp_lrt,
                         exp_time_min=exp_time_tmp)
# } else{
#   stop('Cluster do not match between exponential and spherical models...')
# }

# Merge spatial and non-spatial DE results
res_de = dplyr::full_join(non_sp_de, sp_de, by=c('sample', 'gene', 'cluster1', 'cluster2')) %>%
  dplyr::mutate(spatialmod=dplyr::case_when(is.na(exp_p_val) ~ '2', TRUE ~ '1')) %>% # To put on top of table the genes with spatial models
  dplyr::arrange(spatialmod, cluster1) %>%
  dplyr::select(-spatialmod) %>% 
  tibble::add_column(st_tech=stringr::str_extract(stlist_fp, 'geomx|visium|smi'), .before=1)

rm(non_sp_de, sp_de) # Clean environment

outfp = gsub('data', 'results', stlist_fp)
outfp = gsub('stlist_w_clusters', 'stdiff_results', outfp)
outfp = gsub('.RDS', paste0('_', samplename, '_', genename, '_', cluster1, '_', cluster2, '.csv'), outfp)
write.csv(res_de, file=outfp, row.names=F, quote=F)

