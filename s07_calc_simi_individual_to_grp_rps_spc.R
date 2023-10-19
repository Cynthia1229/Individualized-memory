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

subj_behav <- read_csv(
  path(bids_dir, "facename_behav.csv"), show_col_types = F
) %>% 
  filter(lag == "8_15")

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

grp_avg_rps_spc <- grp_rps_spc %>% 
  semi_join(subj_behav, by = "subj_id") %>% 
  group_by(roi_name, i1_trial_id, i2_trial_id, pair_index) %>% 
  summarise(avg_simi = mean(similarity, na.rm = T), .groups = "drop")

corr_subj_grp_rps_spc <- grp_rps_spc %>% 
  full_join(
    grp_avg_rps_spc, 
    by = c("roi_name", "i1_trial_id", "i2_trial_id", "pair_index")
  ) %>% 
  group_nest(subj_id, roi_name) %>% 
  mutate(
    mdl = map(
      data,
      ~ cor.test(.x$similarity, .x$avg_simi, method = "pearson") %>% 
        broom::tidy()
    )
  ) %>% 
  unnest(mdl) %>% 
  dplyr::select(-data) %>% 
  mutate(fisher_zr = psych::fisherz(estimate))

save(
  grp_avg_rps_spc, corr_subj_grp_rps_spc, 
  file = path(bids_dir, "code", "01_stats_rps_spc", "results", "res_corr_subj_grp_rps_spc.RData")
)


