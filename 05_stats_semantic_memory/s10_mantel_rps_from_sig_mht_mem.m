%%% only calculating Pearson correlation
base_dir = '/brain/guixue/Jintao/eeg_word_memory';
deriv_dir = fullfile(base_dir, 'derivatives');
res_dir = fullfile(base_dir, 'code', '02_mantel_mht_rps_spc', 'results');

% add Mantal_test function
addpath(fullfile(base_dir, 'code', 'code_parameters_sisi', 'Mantel_test_new'))

demean_id = 'nodemean';

% Manhattan distance matrix of memory state
mem_mht = load(fullfile(base_dir, 'word_mht_dist_subj206.mat'));

load(fullfile(deriv_dir, ['rps_spc_' demean_id{1}], 'group', ...
    ['dist_cross_subj_rps_spc_pairwise_' demean_id{1} '.mat']), ...
    'dist_cross_subj_rps_spc'); % region * time * subj * subj
% region 5 and time 600 ~ 100 ms
dist_cross_subj_rps_spc = squeeze(dist_cross_subj_rps_spc(5, 31:end, :,:));
dist_cross_subj_rps_spc = squeeze(nanmean(dist_cross_subj_rps_spc, 1));

