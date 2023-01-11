%%% generate residual matrix

bids_dir = 'D:\Workplace\shared_memory';
data_dir = fullfile(bids_dir, 'code', '05_stats_semantic_memory');
% deriv_dir = fullfile(bids_dir, 'derivatives');
% behav_dir = fullfile(bids_dir, 'code', '00_behav_analysis');
res_dir = fullfile(data_dir, 'results');

% n_perms = 10;

demean_id = 'nodemean';

% add Mantal_test function
% addpath(fullfile(bids_dir, 'code', 'JPMantel_partial'))

subj_behav = readtable(fullfile(data_dir, 'subj_word_behav_infos_206.csv'));

subj_ctrl = subj_behav(subj_behav.hitrate_rk <= 1 & subj_behav.hitrate_rk >= 0, :);
n_subjs = height(subj_ctrl);
idx_tril = find(tril(ones(n_subjs), -1));

idx_subj = ismember(subj_behav.subj_id, subj_ctrl.subj_id);

% behavioral Manhattan distance
confounds = load(fullfile(data_dir, 'data', ...
    ['dat_word_confs_dprime_mat_subj' num2str(n_subjs) '.mat'])); % confounds
confs_list = {'dprime_rk_dist'}; % 'age_dist', 'sex_dist', 

mems_mht = load(fullfile(data_dir, 'data', 'dat_word_mht_dist_subj206.mat'));
mems_list = {'coarse'}; %'precise', 

for imem = 1:length(mems_list)
    mem_id = mems_list{imem};
    mems_mht.(mem_id) = mems_mht.(mem_id)(idx_subj, idx_subj);
end

% load representational space distance
load(fullfile(data_dir, 'data', ...
    ['dist_cross_subj_rps_spc_pairwise_' demean_id '.mat']), ...
    'dist_cross_subj_rps_spc') % region * time * subj * subj

rps_dist_tmp = squeeze(dist_cross_subj_rps_spc(1, 1, :, :));
rps_dist_tmp = rps_dist_tmp(idx_subj, idx_subj);
idx_mat = find(~isnan(rps_dist_tmp));
idx_use = intersect(idx_mat, idx_tril);
%clear rps_dist_tmp

n_regions = size(dist_cross_subj_rps_spc, 1);
n_times = size(dist_cross_subj_rps_spc, 2);

for imem = 1:length(mems_list)
    tic
    memory_id = mems_list{imem};
    %fprintf('Mantel test: %s, %s ... \n', demean_id, memory_id)
    mem_mat = mems_mht.(memory_id);
    
    % standardized matrix
    mem_zscore = zscore(mem_mat(idx_use));
    
    for iconfs = 1:length(confs_list)
        
        conf_id = confs_list{iconfs};
        conf_mat = confounds.(conf_id);
        
        % standardized matrix
        conf_zscore = zscore(conf_mat(idx_use));
        
        % linear model
        lm = fitlm(conf_zscore, mem_zscore); % fitlm(X, y)
        
        resid_tmp = zeros(n_subjs);
        resid_tmp(idx_use) = lm.Residuals.Raw;
        
        resid_tmp = resid_tmp + resid_tmp'; % to symmetrical matrix
        resid_tmp(isnan(rps_dist_tmp)) = nan;
        resid_mem.(conf_id) = resid_tmp;
        clear resid_tmp lm
    end
    save(fullfile(data_dir, 'data', ['resid_dprime_' memory_id '_subj' num2str(n_subjs) '.mat']), ...
        'resid_mem')
end

for iconfs = 1:length(confs_list)
    
    conf_id = confs_list{iconfs};
    conf_mat = confounds.(conf_id);
    
    % standardized matrix
    conf_zscore = zscore(conf_mat(idx_use));
    
    %resid_rps_mean_perm = zeros(n_regions, n_times, n_subjs, n_subjs);
    for iregion = 1:n_regions
        for itime = 1:n_times
            
            rps_dist_tmp = squeeze(dist_cross_subj_rps_spc(iregion, itime, :, :));
            rps_dist_tmp = single(rps_dist_tmp(idx_subj, idx_subj));
            
            % standardized matrix
            rps_dist_zscore = zscore(rps_dist_tmp(idx_use));
            % linear model
            lm = fitlm(conf_zscore, rps_dist_zscore); % fitlm(X, y)
            
            resid_tmp = zeros(n_subjs);
            resid_tmp(idx_use) = lm.Residuals.Raw;
            
            resid_tmp = resid_tmp + resid_tmp'; % to symmetrical matrix
            resid_tmp(isnan(rps_dist_tmp)) = nan;
            resid_isps_dist.(conf_id)(iregion, itime, :, :) = single(resid_tmp);
            clear resid_tmp lm
        end
    end
end

save(fullfile(data_dir, 'data', ...
    ['resid_dprime_dist_cross_subj' num2str(n_subjs) '_isps_pairwise_' ...
    demean_id '.mat']), 'resid_isps_dist')

