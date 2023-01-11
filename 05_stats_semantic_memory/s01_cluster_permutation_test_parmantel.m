function s01_cluster_permutation_test_parmantel(shuffle_num)

rng(shuffle_num) % start from selected seed
%clear,clc;
bids_dir = 'D:\Workplace\shared_memory\code\05_stats_semantic_memory';
%     deriv_dir = fullfile(bids_dir, 'derivatives');
% behav_dir = fullfile(bids_dir, 'data');
res_dir = fullfile(bids_dir, 'results');

demean_id = 'nodemean';
confs_list = {'dprime_rk_dist'};

mems_list = {'coarse'};

n_subjs = 186;

% add Mantal test function
addpath(fullfile(bids_dir, 'Mantel_test_new'))

% load representational space distance
load(fullfile(bids_dir, 'data', ...
    ['resid_dprime_dist_cross_subj' num2str(n_subjs) '_isps_pairwise_' ...
    demean_id '.mat']), 'resid_isps_dist') % structure
% region * time * subj * subj
n_regions = size(resid_isps_dist.(confs_list{1}), 1);
n_times = size(resid_isps_dist.(confs_list{1}), 2);

for imem = 1:length(mems_list)
    memory_id = mems_list{imem};
    
    % get output name
    [~, out_name, ~] = fileparts(fullfile(res_dir, ...
        ['res_parmantel_dprime_mht-word-' memory_id '_isps-word_subjs' ...
        num2str(n_subjs) '_' demean_id '.mat']));
    
    load(fullfile(bids_dir, 'data', ['resid_dprime_' memory_id '_subj' num2str(n_subjs) ...
        '.mat']), 'resid_mem')
    
    %%%%% RUN Permutation %%%%
    
    % shuffle_num = 1000;
    % shuffle_num = 1;
    
    for iconfs = 1:length(confs_list)
        conf_id = confs_list{iconfs};
        
        mem_mat = resid_mem.(conf_id);
        neural_mat = resid_isps_dist.(conf_id);
        
        for ishuffle = 1:shuffle_num
            
            tic
            fprintf('Mantel test: %s, %s, shuffle%d ... \n', demean_id, memory_id, ishuffle)
            % shuffle memory state distance (Manhattan distance)
            ipe = randperm(size(mem_mat, 1));
            mht_temp = mem_mat(ipe, ipe);
            
            % Calculate rho values for correlation
            Corr_beh_STPS_MHT_tem_rho = zeros(n_regions, n_times); % 6 = number of regions
            Corr_beh_STPS_MHT_tem_MTpval = zeros(n_regions, n_times); % 6 = number of regions
            
            for iregion = 1:n_regions % only permute original significant regions
                for itime = 1:n_times
                    mtx1 = squeeze(neural_mat(iregion,itime,:,:));
                    [mrho, mpval] = bramila_mantel(mtx1, mht_temp, 1000, 'pearson');
                    Corr_beh_STPS_MHT_tem_rho(iregion, itime) = mrho;
                    Corr_beh_STPS_MHT_tem_MTpval(iregion, itime) = mpval;
                end
            end
            
            % Get the significant cluster
            corr_Hmap_tem = Corr_beh_STPS_MHT_tem_MTpval <= 0.05;
            
            % Extract significant clusters
            rsum_max_tem_MHT = zeros(n_regions, 1); % 6 = number of regions
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
                        r_sigclu_tem_sum(ic) = sum(Corr_beh_STPS_MHT_tem_rho(iregion, iloc));
                    end
                    % extract biggest cluster for each region
                    rsum_max_tem_MHT(iregion, 1) = max(r_sigclu_tem_sum);
                end
                clear ic L cluster_value
            end
            
            % extract biggest cluster across all regions
            rsum_max_mht.(conf_id)(ishuffle) = max(rsum_max_tem_MHT);
        end
        
        out_dir = fullfile(res_dir, ['parmantel_dprime_perm_' memory_id]);
        
        if ~ exist(out_dir, 'dir')
            mkdir(out_dir)
        end
    end
    save(fullfile(out_dir, [out_name '_rsum_max_all' num2str(shuffle_num, '%04d') '.mat']), ...
        'rsum_max_mht')
    toc
end


