bids_dir = '/seastor/jintaosheng1229/eeg_word_memory';

dat_dir = fullfile(bids_dir, 'derivatives', 'rps_spc_nodemean', 'group');
demean_id = 'nodemean';

load(fullfile(dat_dir, 'group_rps_spc_nodemean.mat'))
mean_grp_rps = squeeze(nanmean(grp_rps_spc, 4));
for isub = 1:size(grp_rps_spc, 4)
    sub_rps_spc = grp_rps_spc(:,:,:,isub);
    for i = 1:size(grp_rps_spc, 1)
        for j = 1:size(grp_rps_spc, 2)
            R(i,j, isub) = corr(squeeze(sub_rps_spc(i,j,:)), ...
                squeeze(mean_grp_rps(i,j,:)), 'rows', 'pairwise');
        end
    end
end
save(fullfile(bids_dir, 'code', '01_model_rsa', 'results', ...
    ['corr_grp_sub_rps_spc_' demean_id '.mat']), 'R')

%% average RPS and hitrate
clear; clc
bids_dir = '/seastor/jintaosheng1229/eeg_word_memory';
behav_subjs = readtable(fullfile(bids_dir, 'subj_word_behav_infos_206.csv'));

load(fullfile(bids_dir, 'code', '01_model_rsa', 'results', ...
    ['corr_grp_sub_rps_spc_' demean_id '.mat']))

fisherz = atanh(R);

for i = 1:size(grp_rps_spc, 1)
    for j = 1:size(grp_rps_spc, 2)
        [R_real(i,j), P_real(i,j)] = corr(behav_subjs.hitrate_rk, squeeze(fisherz(i,j,:)));
    end
end


n_regions = size(R, 1);
n_times = size(R, 2);

shuffle_num = 1000;
%shuffle_num = 1;

hitrate = behav_subjs.hitrate_rk;
%mht_mat = getfield(mem_mht, memory_id);

rsum_max_mean = zeros(shuffle_num, 1);

for ishuffle = 1:shuffle_num
    
    % shuffle memory state distance (Manhattan distance)
    ipe = randperm(size(hitrate, 1));
    hr_temp = hitrate(ipe);
    
    % Calculate rho values for correlation
    Corr_beh_STPS_mean_tem_rho = zeros(n_regions, n_times); % 6 = number of regions
    Corr_beh_STPS_mean_tem_MTpval = zeros(n_regions, n_times); % 6 = number of regions
    
    for iregion = 1:n_regions % only permute original significant regions
        for itime = 1:n_times
            [R_behav_rps, P_behav_rps] = corr(squeeze(fisherz(iregion,itime,:)), ...
                hr_temp, 'rows', 'pairwise');
            Corr_beh_STPS_mean_tem_rho(iregion, itime) = R_behav_rps;
            Corr_beh_STPS_mean_tem_MTpval(iregion, itime) = P_behav_rps;
        end
    end
    
    % Get the significant cluster
    corr_Hmap_tem = Corr_beh_STPS_mean_tem_MTpval <= 0.05;
    
    % Extract significant clusters
    rsum_max_tem_mean = zeros(n_regions, 1); % 6 = number of regions
    for iregion = 1:n_regions
        [L, ~] = bwlabeln(corr_Hmap_tem(iregion,:), 8);
        cluster_value = unique(L);
        cluster_value = setdiff(cluster_value, 0); % remove o
        num_clusters = length(cluster_value);
        if ~isempty(cluster_value)
            % calculate the original r-sum for significant clusters
            r_sigclu_tem_sum = zeros(num_clusters, 1);
            for ic = 1:num_clusters
                iloc = L == cluster_value(ic);
                r_sigclu_tem_sum(ic) = sum(Corr_beh_STPS_mean_tem_rho(iregion, iloc));
            end
            % extract biggest cluster for each region
            rsum_max_tem_mean(iregion, 1) = max(r_sigclu_tem_sum);
        end
        clear ic L cluster_value
    end
    
    % extract biggest cluster across all regions
    rsum_max_mean(ishuffle) = max(rsum_max_tem_mean);
end

save(fullfile(bids_dir, 'code', '01_model_rsa', 'results', ...
    ['corr_grp_sub_rps_spc_hitrate_subjs206_' demean_id '.mat']), ...
    'R_real', 'P_real', 'rsum_max_mean')

