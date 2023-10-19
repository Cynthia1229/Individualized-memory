# calculating representational space similarity across participants


# Prepare----
# Load packages
pacman::p_load(tidyverse, here, fs, glue)

# Directories
bids_dir <- path(here())
deriv_dir <- path(bids_dir, "derivatives")

desc_id <- "leastsquare-separate-hp0p01-unsmooth-tstats"
space_id <- "T1w-to-MNI"
demean_id <- "nodemean"

data_dir <- path(
  deriv_dir, "representational_space", 
  str_glue("rois_{demean_id}_all-wi-bi_thres3p1_unsmo_{demean_id}")
)

message("Loading representational space ... ")

load(
  file = path(
    data_dir, "group", 
    str_glue("group_space-{space_id}_desc-{desc_id}-zr_similarity.RData")
  )
)

subj_rs_epi <- read.table(path(bids_dir, "subjs_fn_rest_461.txt")) %>%
  rename(subj_id = V1)

subj_behav <- read_csv(
  path(bids_dir, "facename_behav.csv"), show_col_types = F
  ) %>% 
  filter(lag == "8_15") %>% 
  #filter(FNRecog_hit <= 0.75, FNRecog_hit >= 0.25) 
  semi_join(subj_rs_epi, by = "subj_id")
  
n_subj <- nrow(subj_behav)

message("Calculating similarity ... ")

simi_rps_spc_cross_sub_norep <- grp_rps_spc %>% 
  semi_join(subj_behav, by = "subj_id") %>% 
  select(subj_id, roi_name, similarity, pair_index) %>% 
  pivot_wider(names_from = "subj_id", values_from = "similarity") %>% 
  select(-pair_index) %>% 
  group_nest(roi_name) %>% 
  mutate(
    res_corr = map(
      data, ~ cor(.x, method = "pearson", use = "complete.obs")
    )
  ) %>% 
  select(-data) 

message("Saving out ... ")
# save outfile
save(
  simi_rps_spc_cross_sub_norep, 
  file = path(
    data_dir, "group", 
    str_glue("simi_rps_spc_cross_sub_norep_subj{n_subj}_{demean_id}.RData")
  )
)

