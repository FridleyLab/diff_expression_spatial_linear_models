##
# Prepare STlist for STdiff analysis in HPC (GeoMx - Liver Spatial Organ Atlas)

library('tidyverse')
library('spatialGE')

# Specify spreadsheet file path
# Download from: https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/human-liver/
liver_fp = '../data/geomx_spatial_organ_atlas_data/hu_liver_count_results/Export3_BiologicalProbeQC.xlsx'

# Read spreadsheet with ROI meta data and clean ROI names to match with
# expression data after using janitor::clean_names
segm_df = readxl::read_excel(liver_fp, sheet=1) %>%
  janitor::clean_names() %>%
  filter(slide_name %in% c("hu_liver_001", "hu_liver_002")) %>%
  mutate(segment_display_name=tolower(segment_display_name) %>%
           str_replace_all(., "[ \\-\\|]+", "_") %>%
           str_replace_all(., "_+", "_"))

# Read expression data and save in individual data frames for each slide
expr_df1 = readxl::read_excel(liver_fp, sheet=3) %>%
  janitor::clean_names() %>%
  column_to_rownames('target_name') %>%
  select(segm_df[['segment_display_name']][segm_df[['slide_name']] == "hu_liver_001"]) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('segment_display_name') %>%
  mutate(segment_display_name=str_extract(segment_display_name, "mipcc_dsp_a_[0-9]+c_01_[0-9]{3}")) %>% # Remove segment name
  group_by(segment_display_name) %>%
  summarise_all(sum) %>%
  column_to_rownames('segment_display_name') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('genename')

expr_df2 = readxl::read_excel(liver_fp, sheet=3) %>%
  janitor::clean_names() %>%
  column_to_rownames('target_name') %>%
  select(segm_df[['segment_display_name']][segm_df[['slide_name']] == "hu_liver_002"]) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('segment_display_name') %>%
  mutate(segment_display_name=str_extract(segment_display_name, "mipcc_dsp_a_[0-9]+c_01_[0-9]{3}")) %>%
  group_by(segment_display_name) %>%
  summarise_all(sum) %>%
  column_to_rownames('segment_display_name') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('genename')

rm(liver_fp) # Clean env

# Prepare coordinate data
coor_df1 = segm_df %>%
  filter(slide_name == "hu_liver_001") %>%
  select(segment_display_name, roi_coordinate_y, roi_coordinate_x) %>%
  mutate(segment_display_name=str_extract(segment_display_name, "mipcc_dsp_a_[0-9]+c_01_[0-9]{3}")) %>%
  distinct()

coor_df2 = segm_df %>%
  filter(slide_name == "hu_liver_002") %>%
  select(segment_display_name, roi_coordinate_y, roi_coordinate_x) %>%
  mutate(segment_display_name=str_extract(segment_display_name, "mipcc_dsp_a_[0-9]+c_01_[0-9]{3}")) %>%
  distinct()

# Create STlist
geomx = STlist(rnacounts=list(hu_liver_001=expr_df1, hu_liver_002=expr_df2),
               spotcoords=list(hu_liver_001=coor_df1, hu_liver_002=coor_df2))
geomx = transform_data(geomx)

# Add layer annotations to data sets
for(i in c('hu_liver_001', 'hu_liver_002')){
  df_tmp = segm_df %>%
    filter(slide_name == i) %>%
    mutate(segment_display_name=str_extract(segment_display_name, "mipcc_dsp_a_[0-9]+c_01_[0-9]{3}")) %>%
    select(segment_display_name, annots=type) %>%
    mutate(annots=tolower(annots) %>%
             str_replace_all(., '[ \\/]+', "_") %>%
             str_replace_all(., "_+", "_")) %>%
    mutate(annots=case_when(str_detect(annots, "cv_") ~ 'central_vein',
                            str_detect(annots, "pv_") ~ 'portal_vein',
                            str_detect(annots, "liver_") ~ 'liver_sin_mac',
                            TRUE ~ annots)) %>%
    distinct()

  geomx@spatial_meta[[i]] = geomx@spatial_meta[[i]] %>%
    left_join(., df_tmp, by=c('libname'='segment_display_name'))

  rm(df_tmp) # clean env
}

rm(segm_df) # Clean env

# Extract genes along a gradient of average expression\
# Make combinations of gene x cluster
tmp_df = list()
for(i in c('hu_liver_001', 'hu_liver_002')){
  genes = geomx@gene_meta[[i]] %>% arrange(desc(gene_mean))
  set.seed(12345)
  genes = sample(genes[['gene']], 5000)
  meta = unique(geomx@spatial_meta[[i]][['annots']])
  tmp_df[[i]] = expand_grid(genes, meta, '../data/geomx_stlist_w_clusters_soaliver.RDS', i)
  rm(meta, genes) # Clean env
}
tmp_df = do.call(bind_rows, tmp_df)
write.table(tmp_df, col.names=F, row.names=F, quote=F,
            file='diff_expr_hpcrun_spamm/data/geomx_gene_meta_combo_soaliver.txt')

rm(tmp_df) # Clean env

saveRDS(geomx, 'diff_expr_hpcrun_spamm/data/geomx_stlist_w_clusters_soaliver.RDS')

