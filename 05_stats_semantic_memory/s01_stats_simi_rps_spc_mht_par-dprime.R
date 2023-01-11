# Mantel test

# Prepare----
# Load packages
pacman::p_load(tidyverse, here, fs, glue, R.matlab)

bids_dir <- path("D:/Workplace/shared_memory/code/05_stats_semantic_memory")
# deriv_dir <- path(bids_dir, "derivatives")

demean_id <- "nodemean"

# rsa_dir <- path(deriv_dir, str_glue("rps_spc_{demean_id}"), "group")

# Manhattan distance
load(
  file = path(
    bids_dir, "data", "dat_mht_distance_word_subj206.RData"
  )
)
mht_dist_mat <- mht_dist_mat %>% mutate(subj = 206)

# Confounds distance
load(
  file = path(
    bids_dir, "data", "dat_confs_dprime_mat_subj206.RData"
  )
)
confounds_pairs_mat <- confounds_pairs_mat %>% mutate(subj = 206)

# Representational space distance
neural4d <- readMat(
  path(
    bids_dir, "data", str_glue("dist_cross_subj_rps_spc_pairwise_{demean_id}.mat")
  )
) %>% 
  pluck("dist.cross.subj.rps.spc")

neural_dfs = tibble() # data.frame() is not good
for (iregion in seq(1,6)) {
  for (itime in seq(1,47)) {
    neural2d <- neural4d[iregion,itime,,]
    neural_df <- tibble(region = iregion, time = itime, neural_mat = list(neural2d))
    neural_dfs <- bind_rows(neural_dfs, neural_df)
  }
}

# partial mantel test
res_parmantel_mht_rps_spc_coarse <- neural_dfs %>% 
  #filter(region == 1, time == 1) %>% 
  mutate(subj = 206) %>%
  full_join(confounds_pairs_mat, by = "subj") %>% 
  mutate(
    mdl = map2(
      neural_mat, mat_dist, 
      ~ vegan::mantel.partial(
        .x, mht_dist_mat$mht_mat[[1]], .y, method = "pearson", 
        permutations = 9999, na.rm = TRUE
      )
    ), 
    statistic_r = map_dbl(mdl, ~ pluck(.x, "statistic")),
    p.value = map_dbl(mdl, ~ pluck(.x, "signif"))
  )

save(
  res_parmantel_mht_rps_spc_coarse,
  file = path(
    bids_dir, "results", 
    str_glue("res_parmantel-dprime_mht-word-coarse_rps-spc-word_subj206_{demean_id}.RData")
  )
)

