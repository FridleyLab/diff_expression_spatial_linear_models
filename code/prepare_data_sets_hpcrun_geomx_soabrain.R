##
# Prepare STlist for STdiff analysis in HPC (GeoMx - Brain cortex Spatial Organ Atlas)

library('tidyverse')
library('spatialGE')

# Specify spreadsheet file path
# Download from: https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/human-brain/
brain_fp = '../data/geomx_spatial_organ_atlas_data/hu_brain_count_results/Export3_BiologicalProbeQC.xlsx'

# Read spreadsheet with ROI meta data and clean ROI names to match with expression data
# Select only cortex ROIs
segm_df = readxl::read_excel(brain_fp, sheet=1) %>%
  janitor::clean_names() %>%
  filter(slide_name %in% c("hu_brain_001", "hu_brain_004a")) %>%
  filter(region == "Cortex") %>%
  mutate(segment_display_name=tolower(segment_display_name) %>%
           str_replace_all(., "[ \\|]+", "_") %>%
           str_replace_all(., "_+", "_") %>%
           str_replace(., "^7", "x7"))

# Read expression data and save in individual data frames for each slide
expr_df1 = readxl::read_excel(brain_fp, sheet=3) %>%
  janitor::clean_names() %>%
  column_to_rownames('target_name') %>%
  select(segm_df[['segment_display_name']][segm_df[['slide_name']] == "hu_brain_001"]) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('segment_display_name') %>%
  mutate(segment_display_name=str_extract(segment_display_name, "x74_[0-9]{3}")) %>% # Remove segment name
  group_by(segment_display_name) %>%
  summarise_all(sum) %>%
  column_to_rownames('segment_display_name') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('genename') %>%
  filter(!str_detect(genename, 'Neg')) # Remove negative probe counts

expr_df2 = readxl::read_excel(brain_fp, sheet=3) %>%
  janitor::clean_names() %>%
  column_to_rownames('target_name') %>%
  select(segm_df[['segment_display_name']][segm_df[['slide_name']] == "hu_brain_004a"]) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('segment_display_name') %>%
  mutate(segment_display_name=str_extract(segment_display_name, "cortex_29_[0-9]{3}")) %>% # Remove segment name
  group_by(segment_display_name) %>%
  summarise_all(sum) %>%
  column_to_rownames('segment_display_name') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('genename') %>%
  filter(!str_detect(genename, 'Neg')) # Remove negative probe counts

rm(brain_fp) # Clean env

# Prepare coordinate data
coor_df1 = segm_df %>%
  filter(slide_name == "hu_brain_001") %>%
  select(segment_display_name, roi_coordinate_y, roi_coordinate_x) %>%
  mutate(segment_display_name=str_extract(segment_display_name, "x74_[0-9]{3}")) %>%
  distinct()

coor_df2 = segm_df %>%
  filter(slide_name == "hu_brain_004a") %>%
  select(segment_display_name, roi_coordinate_y, roi_coordinate_x) %>%
  mutate(segment_display_name=str_extract(segment_display_name, "cortex_29_[0-9]{3}")) %>%
  distinct()

# Create STlist
geomx = STlist(rnacounts=list(hu_brain_001=expr_df1, hu_brain_004a=expr_df2),
               spotcoords=list(hu_brain_001=coor_df1, hu_brain_004a=coor_df2))
geomx = transform_data(geomx)

# Add layer annotations to data sets
for(i in c('hu_brain_001', 'hu_brain_004a')){
  df_tmp = segm_df %>%
    filter(slide_name == i) %>%
    mutate(segment_display_name=str_extract(segment_display_name, "x74_[0-9]{3}|cortex_29_[0-9]{3}")) %>%
    select(segment_display_name, layer=group) %>%
    mutate(annots=tolower(layer) %>%
             str_replace_all(., '[ \\/]+', "_") %>%
             str_replace_all(., "_+", "_")) %>%
    distinct()

  geomx@spatial_meta[[i]] = geomx@spatial_meta[[i]] %>%
    left_join(., df_tmp, by=c('libname'='segment_display_name'))

  rm(df_tmp) # clean env
}

rm(segm_df) # Clean env

# Extract genes along a gradient of average expression\
# Make combinations of gene x cluster
tmp_df = list()
for(i in c('hu_brain_001', 'hu_brain_004a')){
  genes = geomx@gene_meta[[i]] %>% arrange(desc(gene_mean))
  set.seed(12345)
  genes = sample(genes[['gene']], 5000)
  meta = unique(geomx@spatial_meta[[i]][['annots']])
  tmp_df[[i]] = expand_grid(genes, meta, '../data/geomx_stlist_w_clusters_soabrain.RDS', i)
  rm(meta, genes) # Clean env
}
tmp_df = do.call(bind_rows, tmp_df)
write.table(tmp_df, col.names=F, row.names=F, quote=F,
            file='diff_expr_hpcrun_spamm/data/geomx_gene_meta_combo_soabrain.txt')

rm(tmp_df) # Clean env

saveRDS(geomx, 'diff_expr_hpcrun_spamm/data/geomx_stlist_w_clusters_soabrain.RDS')

