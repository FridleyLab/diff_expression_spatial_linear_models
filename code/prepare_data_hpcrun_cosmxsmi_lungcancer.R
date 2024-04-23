##
# Prepare STlist for STdiff analysis in HPC (CosmxSMI - Non small cell lung cancer)

library('tidyverse')
library('spatialGE')

# Read SMI data before passing to STlist
# Download from: https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/
load('../data/cosmx_smi_lung_nsclc_nanostring/data/Processed Data Giotto Object/SMI_Giotto_Object.RData')

# Create list to store data and pass to STlist
expr_dfs = list()
coor_dfs = list()

# Organize data for slide Lung5-3, FOV 28
# Filter cell types with poor representation (keep if > 10 cells)
cellids = gem@cell_metadata$rna %>%
  filter(Slide_name == 'Lung5-3' & fov == '28') %>%
  mutate(cell_type=tolower(cell_type) %>%
           str_replace_all(., '[ \\-]', '_')) %>%
  filter(!(cell_type %in% c('treg', 't_cd8_memory', 'tumor_13')))
expr_dfs[['lung5rep3_fov28']] = as.data.frame(gem@expression$rna$raw[, cellids[['cell_ID']]]) %>%
  rownames_to_column(var='genename')
coor_dfs[['lung5rep3_fov28']] = gem@spatial_locs$viz[gem@spatial_locs$viz$cell_ID %in% cellids[['cell_ID']], ] %>%
  select(cell_ID, sdimy, sdimx)

# Organize data for slide Lung6, FOV 20
cellids = gem@cell_metadata$rna %>%
  filter(Slide_name == 'Lung6' & fov == '20') %>%
  mutate(cell_type=tolower(cell_type) %>%
           str_replace_all(., '[ \\-]', '_')) %>%
  filter(!(cell_type %in% c( 'b_cell', 't_cd4_naive', 'mast', 'plasmablast', 'tumor_13', 'tumor_5', 'tumor_9')))
expr_dfs[['lung6_fov20']] = as.data.frame(gem@expression$rna$raw[, cellids[['cell_ID']]]) %>%
  rownames_to_column(var='genename')
coor_dfs[['lung6_fov20']] = gem@spatial_locs$viz[gem@spatial_locs$viz$cell_ID %in% cellids[['cell_ID']], ] %>%
  select(cell_ID, sdimy, sdimx)

rm(cellids) # Clear env

# Create STlist from SMI data sets (this step to easily extract count matrices)
smi = STlist(rnacounts=expr_dfs, spotcoords=coor_dfs)
smi = transform_data(smi)

# Add tissue niche annotations to STlist
for(i in c('lung5rep3_fov28', 'lung6_fov20')){
  df_tmp = gem@cell_metadata$rna %>%
    filter(cell_ID %in% smi@spatial_meta[[i]][['libname']]) %>%
    mutate(cell_type=tolower(cell_type) %>%
           str_replace_all(., '[ \\-]', '_')) %>%
    #select(cell_ID, cell_type, niche)
    select(cell_ID, annots=cell_type)

  smi@spatial_meta[[i]] = smi@spatial_meta[[i]] %>%
    left_join(., df_tmp, by=c('libname'='cell_ID'))

  rm(df_tmp) # clean env
}

rm(gem) # Clean env

# Make combinations of gene x cluster
# Get all genes in data set
tmp_df = list()
for(i in c('lung5rep3_fov28', 'lung6_fov20')){
  genes = smi@gene_meta[[i]] %>% arrange(desc(gene_mean))
  genes = genes[['gene']]
  meta = unique(smi@spatial_meta[[i]][['annots']])
  tmp_df[[i]] = expand_grid(genes, meta, '../data/smi_stlist_w_clusters_lungcancer.RDS', i)
  rm(meta, genes) # Clean env
}
tmp_df = do.call(bind_rows, tmp_df)
write.table(tmp_df, col.names=F, row.names=F, quote=F,
            file='diff_expr_hpcrun_spamm/data/smi_gene_meta_combo_lungcancer.txt')

rm(tmp_df) # Clean env

saveRDS(smi, 'diff_expr_hpcrun_spamm/data/smi_stlist_w_clusters_lungcancer.RDS')

