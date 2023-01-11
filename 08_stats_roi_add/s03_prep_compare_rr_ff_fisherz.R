# setup
pacman::p_load(tidyverse, here, fs, ggsci, ggpubr)

# preparing
bids_dir <- path("/seastor/jintaosheng1229/shared_memory")
deriv_dir <- path(bids_dir, "derivatives")
dat_dir <- path(deriv_dir, "within-item_cross-subject_PS", "group")

demean_id <- "nodemean"

roi_names <- c('PMC', 'vmPFC', 'lHIP', 'rHIP')

subj_stats <- read_csv(
  path(bids_dir, "facename_behav.csv"), show_col_types = FALSE
) %>% 
  filter(lag == "8_15", mfd <= 0.3)

n_subjs <- nrow(subj_stats)

grp_stats <- read_csv(
  path(dat_dir, "group_behav_labels_subj415.csv"), show_col_types = F
) %>% 
  semi_join(subj_stats, by = "subj_id") %>% 
  filter(num_pres == 1) %>% 
  select(subj_id, face_id, memory_score)

simi_wics_mean <- tibble()

for (roi_name in roi_names) {
  simi_wics <- read_csv(
    path(dat_dir, str_glue("subjpairs_correlation_within-item_cross-subj_allitems_{roi_name}.csv")),
    show_col_types = F
  ) %>% 
    left_join(grp_stats, by = c("subj_id1" = "subj_id", "face_id")) %>% 
    rename(i1_memory_score = memory_score) %>% 
    left_join(grp_stats, by = c("subj_id2" = "subj_id", "face_id")) %>% 
    rename(i2_memory_score = memory_score) %>% 
    filter(i1_memory_score == i2_memory_score) %>% 
    mutate(similarity = psych::fisherz(similarity)) %>% 
    group_by(subj_id1, subj_id2, memory_score = i1_memory_score) %>% 
    summarise(similarity = mean(similarity, na.rm = T), .groups = "drop") %>% 
    mutate(roi_name = roi_name)
  simi_wics_mean <- rbind(simi_wics_mean, simi_wics)
}

save(
  simi_wics_mean, 
  file = path(
    bids_dir, "code", "08_stats_roi_add", "data", 
    str_glue("simi_wics_mean_subjpairs_roisadd_{demean_id}_fisherz.RData")
  )
)
