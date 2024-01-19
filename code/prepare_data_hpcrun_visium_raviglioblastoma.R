##
# Prepare STlist for STdiff analysis in HPC (Visium - Glioblastoma)

library('tidyverse')
library('SPATA2')
library('spatialGE')

# Download SPATA objects with spot annotations
SPATAData::downloadSpataObjects(sample_names=c('243_T', '275_T'), overwrite=T,
                                files=c('../data/ravi_gbm_spata_objects/UKF243_T_ST.RDS',
                                        '../data/ravi_gbm_spata_objects/UKF275_T_ST.RDS'))

# Create STlist from Visium data sets (this step to easily extract count matrices)
# Download Visium "outs" folders from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.h70rxwdmj
fps = c('../data/ravi_2022_glioblastoma/UKF243_T_ST/outs/',
        '../data/ravi_2022_glioblastoma/UKF275_T_ST/outs/')

samples = c(
  'UKF243',
  'UKF275')

visium = STlist(rnacounts=fps, samples=samples)
visium = transform_data(visium)

rm(fps, samples) # Clean env

# Add spot annotations from SPATA object (Ravi et al. 2022)
annots_ls = list(UKF243=readRDS('../data/ravi_gbm_spata_objects/UKF243_T_ST.RDS'),
                 UKF275=readRDS('../data/ravi_gbm_spata_objects/UKF275_T_ST.RDS'))

for(i in c('UKF243', 'UKF275')){
  df_tmp = annots_ls[[i]]@fdata[[1]]

  visium@spatial_meta[[i]] = visium@spatial_meta[[i]] %>%
    left_join(., df_tmp %>%
                select(barcodes, annots=segmentation), by=c('libname'='barcodes'))

  rm(df_tmp) # clean env
}

rm(annots_ls) # Clean env

# Extract genes aling a gradient of average expression
# Make combinations of gene x cluster
tmp_df = list()
for(i in c('UKF243', 'UKF275')){
  genes = visium@gene_meta[[i]] %>% arrange(desc(gene_mean))
  set.seed(12345)
  genes = sample(genes[['gene']], 5000)
  meta = unique(visium@spatial_meta[[i]][['annots']])
  tmp_df[[i]] = expand_grid(genes, meta, '../data/visium_stlist_w_clusters_raviglioblastoma.RDS', i)
  rm(meta, genes) # Clean env
}
tmp_df = do.call(bind_rows, tmp_df)
write.table(tmp_df, col.names=F, row.names=F, quote=F,
            file='diff_expr_hpcrun_spamm/data/visium_gene_meta_combo_raviglioblastoma.txt')

rm(tmp_df) # Clean env

saveRDS(visium, 'diff_expr_hpcrun_spamm/data/visium_stlist_w_clusters_raviglioblastoma.RDS')

