# setup
pacman::p_load(tidyverse, here, fs, ggsci, ggpubr, vegan)

# preparing
bids_dir <- path(here())
deriv_dir <- path(bids_dir, "derivatives")

desc_id <- "leastsquare-separate-hp0p01-unsmooth-tstats"
space_id <- "T1w-to-MNI"
demean_id <- "nodemean"

data_dir <- path(
  deriv_dir, "representational_space", 
  str_glue("rois_{demean_id}_all-wi-bi_thres3p1_unsmo_{demean_id}")
)

stats_list <- c("subjs415")

for (istats in stats_list) {
  
  # load behavioral data
  # behavs_smc_mat, confounds_pairs_mat, behavs_smc
  load(
    file = path(
      bids_dir, "code", "00_behav_analysis", 
      str_glue("dat_behavs_mat_{istats}.RData")
    )
  )
  # load manhattan distance
  load(
    file = path(
      bids_dir, "code", "00_behav_analysis", 
      str_glue("dat_mht_distance_{istats}.RData")
    )
  )
  
  mht_dist_mat <- mht_dist_mat %>% mutate(merge_id = "all")
  
  confounds_pairs_mat <- confounds_pairs_mat %>% mutate(merge_id = "all")
  
  load(
    file = path(
      data_dir, "group", 
      str_glue("simi_resid_rps_spc_cross_sub_norep_{istats}_{demean_id}.RData")
    )
  )
  
  message(str_glue("{istats}: Mantel test ... "))
  
  # without confounding
  res_mantel_mht_rps_spc <- simi_resid_rps_spc_cross_sub_norep %>%
    mutate(demean_id = demean_id) %>% 
    mutate(merge_id = "all") %>% 
    full_join(mht_dist_mat, by = "merge_id") %>% 
    mutate(
      mdl = map2(
        res_corr, mht_mat, 
        ~ mantel(
          1 - .x, .y, method = "pearson", permutations = 9999, na.rm = TRUE
        )
      ),
      statistic_r = map_dbl(mdl, ~ pluck(.x, "statistic")),
      p.value = map_dbl(mdl, ~ pluck(.x, "signif"))
    )
  
  message(str_glue("{istats}: Saving out Mantel test ... "))
  save(
    res_mantel_mht_rps_spc,
    file = path(
      bids_dir, "code", "01_stats_rps_spc", "results", 
      str_glue("res_mantel_mht-fn_rps-spc-fn-resid_{istats}_{demean_id}.RData")
    )
  )
  
  message(str_glue("{istats}: Partial Mantel test for association memory ... "))
  
  # with confoundings
  ### association memory
  tmp <- confounds_pairs_mat %>% 
    filter(confounds == "FNRecall_hit" | confounds == "sex" | confounds == "age")
  
  res_parmantel_mht_rps_spc_association <- simi_resid_rps_spc_cross_sub_norep %>%
    mutate(demean_id = demean_id) %>% 
    mutate(merge_id = "all") %>% 
    full_join(tmp, by = "merge_id") %>% 
    mutate(
      mdl = map2(
        res_corr, mat_dist, 
        ~ mantel.partial(
          1 - .x, mht_dist_mat$mht_mat[[1]], .y, method = "pearson", 
          permutations = 9999, na.rm = TRUE
        )
      ), 
      statistic_r = map_dbl(mdl, ~ pluck(.x, "statistic")),
      p.value = map_dbl(mdl, ~ pluck(.x, "signif"))
    )
  
  message(str_glue("{istats}: Partial Mantel test for precise memory ... "))
  
  ### Item memory
  tmp <- confounds_pairs_mat %>% 
    filter(confounds == "FNRecog_hit" | confounds == "sex" | confounds == "age")
  
  res_parmantel_mht_rps_spc_precise <- simi_resid_rps_spc_cross_sub_norep %>%
    mutate(demean_id = demean_id) %>% 
    mutate(merge_id = "all") %>% 
    full_join(tmp, by = "merge_id") %>% 
    mutate(
      mdl = map2(
        res_corr, mat_dist, 
        ~ mantel.partial(
          1 - .x, mht_dist_mat$mht_mat[[2]], .y, method = "pearson", 
          permutations = 9999, na.rm = TRUE
        )
      ), 
      statistic_r = map_dbl(mdl, ~ pluck(.x, "statistic")),
      p.value = map_dbl(mdl, ~ pluck(.x, "signif"))
    )
  
  message(str_glue("{istats}: Saving out partial mantel test ... "))
  
  save(
    res_parmantel_mht_rps_spc_association, 
    res_parmantel_mht_rps_spc_precise,
    file = path(
      bids_dir, "code", "01_stats_rps_spc", "results", 
      str_glue("res_parmantel_mht-fn_rps-spc-fn-resid_{istats}_{demean_id}.RData")
    )
  )
}



