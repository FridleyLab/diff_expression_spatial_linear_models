##
# Prepare STlist for STdiff analysis in HPC (Visium - Brain motor cortex)

library('tidyverse')
library('spatialGE')

# Create STlist from Visium data sets (this step to easily extract count matrices)
# Download Visium data from: https://github.com/LieberInstitute/HumanPilot
fps = c('../data/maynard_2021_prefrontal_cortex/data/151507/',
        '../data/maynard_2021_prefrontal_cortex/data/151673/')

samples = c(
  '151507',
  '151673')

visium = STlist(rnacounts=fps, samples=samples)
visium = transform_data(visium)

rm(fps, samples) # Clean env

# Add layer annotations to Maynard data sets
# Download Visium data from: https://github.com/LieberInstitute/HumanPilot
brain_fps = '../data/maynard_2021_prefrontal_cortex/data/annotated_spots/'
brain_fps = list.files(brain_fps, full.names=T)

for(i in c('151507', '151673')){
  fp_tmp = grep(i, brain_fps, value=T)
  df_tmp = tibble()
  for(fp in fp_tmp){
    layer_tmp = str_extract(fp, "_[LWM0-9]+_barcodes") %>% str_replace_all(., "_|barcodes", "")
    df_tmp = bind_rows(df_tmp,
                       tibble(libname=scan(fp, what='c', quiet=T)) %>%
                         add_column(annots=layer_tmp, .after=1))
  }
  visium@spatial_meta[[i]] = visium@spatial_meta[[i]] %>%
    left_join(., df_tmp, by='libname')

  rm(df_tmp, fp_tmp) # clean env
}

# Remove spots with NA as annotation
for(i in names(visium@spatial_meta)){
  na_spots = visium@spatial_meta[[i]][['libname']][ is.na(visium@spatial_meta[[i]][['annots']]) ]
  visium = filter_data(visium, rm_spots=na_spots, samples=i)

  rm(na_spots) # Clean env
}

# Extract genes along a gradient of average expression\
# Make combinations of gene x cluster
tmp_df = list()
for(i in c('151507', '151673')){
  genes = visium@gene_meta[[i]] %>% arrange(desc(gene_mean))
  set.seed(12345)
  genes = sample(genes[['gene']], 5000)
  meta = unique(visium@spatial_meta[[i]][['annots']])
  tmp_df[[i]] = expand_grid(genes, meta, '../data/visium_stlist_w_clusters_maynardbrain.RDS', i)
  rm(meta, genes) # Clean env
}
tmp_df = do.call(bind_rows, tmp_df)
write.table(tmp_df, col.names=F, row.names=F, quote=F,
            file='diff_expr_hpcrun_spamm/data/visium_gene_meta_combo_maynardbrain.txt')

rm(tmp_df) # Clean env

saveRDS(visium, 'diff_expr_hpcrun_spamm/data/visium_stlist_w_clusters_maynardbrain.RDS')

