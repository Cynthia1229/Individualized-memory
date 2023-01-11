# permutation
pacman::p_load(tidyverse, here, fs, glue, dplyr, R.matlab, vegan, corrr)

# Reading matrix
bids_dir <- path(here())
dat_dir <- path(bids_dir, "derivatives", "semantic_memory")
code_dir <- path(bids_dir, "code", "03_network_analysis")

# Behavirol matrix
load(path(code_dir, "behav-sem_smc_acc-pairs_mat206.RData"))

# Neural matrix
neural4d <- readMat(
  path(dat_dir, "crosssub_pattern_corr_dim_6_46_206_206.mat")
) %>% 
  pluck("crosssub.pattern.corr")


# Preparing neural data
neural_dfs = data.frame()
for (iregion in seq(1,6)) {
  for (itime in seq(1,46)) {
    neural2d <- neural4d[iregion,itime,,]
    neural_df <- tibble(region = iregion, time = itime, neural_mat = list(neural2d))
    neural_dfs <- bind_rows(neural_dfs, neural_df)
  }
}

# Partial Mantel Test
# res_parmantel_semantic <- neural_dfs %>% 
#   mutate(
#     res_mantel = map(
#       neural_mat, 
#       ~ mantel.partial(
#         .x, smc_mat, acc_pairs_mat, 
#         method = "pearson", permutations = 999, na.rm = TRUE
#       )
#     ), 
#     statistic.r = map_dbl(res_mantel, ~ pluck(.x, "statistic")),
#     p.value = map_dbl(res_mantel, ~ pluck(.x, "signif"))
#   ) %>% 
#   select(-c(neural_mat, res_mantel))
# 
# r_parmat <- res_parmantel_semantic %>% 
#   select(-p.value) %>% 
#   pivot_wider(names_from = "time", values_from = "statistic.r")
# write_csv(r_parmat, file = path(dat_dir, "parmantel_smc_rem_know_semantic_r.csv"))
# 
# p_parmat <- res_parmantel_semantic %>% 
#   select(-statistic.r) %>% 
#   pivot_wider(names_from = "time", values_from = "p.value")
# write_csv(p_parmat, file = path(dat_dir, "parmantel_smc_rem_know_semantic_p.csv"))
# 
# write_csv(
#   res_parmantel_semantic, 
#   file = path(dat_dir, "parmantel_smc_rem_know_semantic_real.csv")
# )

smc_longmat <- smc_mat %>% as_cordf() %>% shave() %>% 
  stretch() %>% na.omit()

res_parmantel_semantic_rand1000 <- tibble()

for (iperm in seq(1, 1000)) {
  message(str_glue("Permutation: iperm = {iperm}"))
  smc_randmat <- smc_longmat %>% 
    mutate(shuffler = sample(r)) %>% 
    select(-r) %>% 
    pivot_wider(
      names_from = y, values_from = shuffler, values_fill = 0
    ) %>% 
    select(-x) %>% 
    mutate(`2002` = 0, .before = 1L) %>% 
    as.matrix() %>% 
    rbind(0)
  smc_randmat <- smc_randmat + t(smc_randmat) + diag(nrow(smc_randmat))
  
  res_parmantel_semantic_rand <- neural_dfs %>% 
    mutate(
      res_mantel = map(
        neural_mat, 
        ~ mantel.partial(
          .x, smc_randmat, acc_pairs_mat, 
          method = "pearson", permutations = 999, na.rm = TRUE
        )
      ), 
      statistic.r = map_dbl(res_mantel, ~ pluck(.x, "statistic")),
      p.value = map_dbl(res_mantel, ~ pluck(.x, "signif")),
      perm_times = iperm
    ) %>% 
    select(-c(neural_mat, res_mantel))
  res_parmantel_semantic_rand1000 <- rbind(
    res_parmantel_semantic_rand1000, res_parmantel_semantic_rand
    )
}

# Write out
write_csv(
  res_parmantel_semantic_rand1000, 
  file = path(dat_dir, "parmantel_smc_rem_know_semantic_perm1000.csv")
)

