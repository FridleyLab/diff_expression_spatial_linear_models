##
# Prepare STlist for STdiff analysis in HPC (CosmxSMI - Normal and cancerous liver)

library('tidyverse')
library('spatialGE')

# Read SMI data before passing to STlist
# Download from: https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/
liver_seurat = readRDS('../data/cosmx_smi_liver_nanostring/LiverDataReleaseSeurat_noTranscripts_newUMAP.RDS')


# Subset to selected FOVs
cell_keep_fov1 = liver_seurat@meta.data %>%
  filter(Run_Tissue_name == "NormalLiver" & fov == 174) %>%
  rownames()

cell_keep_fov2 = liver_seurat@meta.data %>%
  filter(Run_Tissue_name == "CancerousLiver" & fov == 366) %>%
  rownames()

# Create list to store data and pass to STlist
expr_dfs = list(NormalLiver_174=liver_seurat@assays$RNA@counts[, cell_keep_fov1],
                CancerousLiver_366=liver_seurat@assays$RNA@counts[, cell_keep_fov2])
coor_dfs = list(NormalLiver_174=liver_seurat@meta.data[cell_keep_fov1, ],
                CancerousLiver_366=liver_seurat@meta.data[cell_keep_fov2, ])

# Organize data
# Filter cell types with poor representation (keep if > 10 cells)
coor_dfs[['NormalLiver_174']] = coor_dfs[['NormalLiver_174']]  %>%
  rownames_to_column('cellid') %>%
  mutate(annots=tolower(cellType) %>%
           str_replace_all(., '[\\+\\.]+', '_')) %>%
  filter(!(annots %in% c('erthyroid_cells', 'portal_endothelial_cells', 'nk_like_cells',
                            'central_venous_lsecs', 'hep_1'))) %>%
  select('cellid', "y_slide_mm", "x_slide_mm")
expr_dfs[['NormalLiver_174']] = expr_dfs[['NormalLiver_174']][, coor_dfs[['NormalLiver_174']][['cellid']] ] %>%
  as.matrix() %>% as.data.frame() %>%
  rownames_to_column(var='genename')

coor_dfs[['CancerousLiver_366']] = coor_dfs[['CancerousLiver_366']]  %>%
  rownames_to_column('cellid') %>%
  mutate(annots=tolower(cellType) %>%
           str_replace_all(., '[\\+\\.]+', '_')) %>%
  filter(!(annots %in% c('cholangiocytes', 'erthyroid_cells', 'gamma_delta_t_cells_1',
                            'nk_like_cells', 'portal_endothelial_cells'))) %>%
  select('cellid', "y_slide_mm", "x_slide_mm")
expr_dfs[['CancerousLiver_366']] = expr_dfs[['CancerousLiver_366']][, coor_dfs[['CancerousLiver_366']][['cellid']] ] %>%
  as.matrix() %>% as.data.frame() %>%
  rownames_to_column(var='genename')

rm(cell_keep_fov1, cell_keep_fov2) # Clear env

# Create STlist from SMI data sets (this step to easily extract count matrices)
smi = STlist(rnacounts=expr_dfs, spotcoords=coor_dfs)
smi = transform_data(smi)

# Add tissue niche annotations to STlist
for(i in c('NormalLiver_174', 'CancerousLiver_366')){
  df_tmp = liver_seurat@meta.data %>%
    rownames_to_column('cellid') %>%
    filter(cellid %in% smi@spatial_meta[[i]][['libname']]) %>%
    mutate(cellType=tolower(cellType) %>%
           str_replace_all(., '[\\+\\.]+', '_')) %>%
    #select(cell_ID, cell_type, niche)
    select(cellid, annots=cellType)

  smi@spatial_meta[[i]] = smi@spatial_meta[[i]] %>%
    left_join(., df_tmp, by=c('libname'='cellid'))

  rm(df_tmp) # clean env
}

# Make combinations of gene x cluster
# Get all genes in data set
tmp_df = list()
for(i in c('NormalLiver_174', 'CancerousLiver_366')){
  genes = smi@gene_meta[[i]] %>% arrange(desc(gene_mean))
  genes = genes[['gene']]
  meta = unique(smi@spatial_meta[[i]][['annots']])
  tmp_df[[i]] = expand_grid(genes, meta, '../data/smi_stlist_w_clusters_liver.RDS', i)
  rm(meta, genes) # Clean env
}
tmp_df = do.call(bind_rows, tmp_df)
write.table(tmp_df, col.names=F, row.names=F, quote=F,
            file='diff_expr_hpcrun_spamm/data/smi_gene_meta_combo_liver.txt')

rm(tmp_df) # Clean env

saveRDS(smi, 'diff_expr_hpcrun_spamm/data/smi_stlist_w_clusters_liver.RDS')

