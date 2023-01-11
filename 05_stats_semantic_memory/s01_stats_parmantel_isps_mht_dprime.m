%%% only calculating Pearson correlation
bids_dir = 'D:\Workplace\shared_memory\code\05_stats_semantic_memory';
% deriv_dir = fullfile(bids_dir, 'derivatives');
% behav_dir = fullfile(bids_dir, 'code', '00_behav_analysis');
res_dir = fullfile(bids_dir, 'results');

n_subjs = 206;

demean_id = 'nodemean';

% add Mantal_test function
addpath(fullfile(bids_dir, 'Mantel_test_new'))

% behavioral Manhattan distance
confs_list = {'dprime_rk_dist'};

mems_list = {'coarse'}; % 'precise', 

% load representational space distance
load(fullfile(bids_dir, 'data', ...
    ['resid_dprime_dist_cross_subj' num2str(n_subjs) '_isps_pairwise_' ...
    demean_id '.mat']), 'resid_isps_dist') % structure
% region * time * subj * subj
n_regions = size(resid_isps_dist.(confs_list{1}), 1);
n_times = size(resid_isps_dist.(confs_list{1}), 2);
%n_subjs = size(resid_rps_mean_perm.age_dist, 3);

for imem = 1:length(mems_list)
    tic
    memory_id = mems_list{imem};
    
    load(fullfile(bids_dir, 'data', ...
        ['resid_dprime_' memory_id '_subj' num2str(n_subjs) '.mat']), ...
        'resid_mem')
    
    fprintf('Mantel test: %s, %s ... \n', demean_id, memory_id)
    
    for iconfs = 1:length(confs_list)
        conf_id = confs_list{iconfs};
        mem_mat = resid_mem.(conf_id);
        neural_mat = resid_isps_dist.(conf_id);
        %pearson_r.(memory_id) = zeros(n_regions, n_times);
        %pvalue.(memory_id) = zeros(n_regions, n_times);
        for iregion = 1:n_regions
            for itime = 1:n_times
                rps_dist_tmp = squeeze(neural_mat(iregion, itime, :, :));
                [pmrval, pmpval] = bramila_mantel(rps_dist_tmp, mem_mat, 10000, 'pearson');
                pearson_r.(conf_id)(iregion, itime) = pmrval;
                pvalue.(conf_id)(iregion, itime) = pmpval;
            end
        end
    end
    
    save(fullfile(res_dir, ['res_parmantel_dprime_mht-word-' memory_id ...
        '_isps-word_subjs' num2str(n_subjs) '_' demean_id '1w.mat']), ...
        'pearson_r', 'pvalue')
    toc
end

