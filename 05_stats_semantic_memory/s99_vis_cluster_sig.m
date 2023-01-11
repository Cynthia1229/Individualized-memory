
%% windows set
n_times = 307; % 1200 ms
n_times = 256; % 1000 ms
window_length = 26; % 100 ms
step = 5; % 20 ms
window_sets = 1:step:(n_times-window_length+1);

141*(1200/307); % 29 time window (141)
211*(1200/307); % 43 time window (211)

%% cluster sig memory and representational space (random 10 times)
clear,clc;
bids_dir = '/brain/guixue/Jintao/eeg_word_memory';
res_dir = fullfile(bids_dir, 'code', '02_mantel_mht_rps_spc', 'results');

addpath(fullfile(bids_dir, 'code')) % add circle_cluster.m function

demean_id = 'nodemean';
n_regions = 6;

memory_types = {'precise', 'coarse'}; % 'precise', 'coarse'

load(fullfile(res_dir, ...
    ['res_mantel_mht-word_rps-spc-word_perm10mean_subjs206_' demean_id '.mat']), ...
    'pearson_r', 'pvalue');

[~, out_name, ~] = fileparts(fullfile(res_dir, ...
    ['res_mantel_mht-word_rps-spc-word_perm10mean_subjs206_' demean_id '.mat']));

for iconf = 1:length(memory_types)
    memory_id = memory_types{iconf};
    
    res_r = pearson_r.(memory_id);
    res_p = pvalue.(memory_id);
    
    % get the original cluster
    corr_Hmap = res_p <= 0.05;
    
    % extract significant clusters
    for iregion = 1:n_regions
        [L, ~] = bwlabeln(corr_Hmap(iregion,:), 8);
        cluster_value = unique(L);
        cluster_value = setdiff(cluster_value, 0);
        if isempty(cluster_value)   % No significant clusters
            Corr_BehNeu_Orig.(memory_id)(iregion).clusters_idx = [];
            Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu = [];
        else
            % calculate the original r-sum for significant clusters
            for ic = 1:length(cluster_value)
                [i] = find(L == cluster_value(ic));
                r_orig_sum_tmp = res_r(iregion,i);
                % get x-axis & y-axis of the significant areas
                Corr_BehNeu_Orig.(memory_id)(iregion).clusters_idx{ic,1} = i;
                Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu(ic) = sum(r_orig_sum_tmp);
                clear i r_orig_sum_tmp
            end
        end
        clear ic L num cluster_value
    end
    
    
    % Compare original_sig_clusters with the perm_max_cluster to calculate p-values
    
    
    % load r sum max
    
    perm_files = dir(fullfile(res_dir, ['perm10mean_' memory_id], 'res_*.mat'));
    n_perms = length(perm_files);
    rmax_perm_all = zeros(n_perms, 1);
    for iperm = 1:n_perms
        load(fullfile(perm_files(iperm).folder, perm_files(iperm).name), 'rsum_max_mht')
        rmax_perm_all(iperm) = rsum_max_mht;
    end
    
    % R cluster
    for iregion = 1:n_regions
        if ~isempty(Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu)
            for ic = 1:length(Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu)
                Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster(ic) = length(find(rmax_perm_all > Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu(ic)))./n_perms;
            end
        else
            Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster = [];
            Corr_BehNeu_Orig.(memory_id)(iregion).p_min_cluster = [];
        end
        
        Corr_BehNeu_Orig.(memory_id)(iregion).p_max_idx = find(Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster <= 0.05);
        Corr_BehNeu_Orig.(memory_id)(iregion).P_cluster_idx = [Corr_BehNeu_Orig.(memory_id)(iregion).p_max_idx];
    end
    
end
save(fullfile(res_dir, ['cluster_' out_name '.mat']), 'Corr_BehNeu_Orig')


%% cluster sig memory and representational space (206 raw data)
clear,clc;
bids_dir = 'D:\Workplace\shared_memory';
res_dir = fullfile(bids_dir, 'code', '05_stats_semantic_memory', 'results');

addpath(fullfile(bids_dir, 'code', '05_stats_semantic_memory')) % add circle_cluster.m function

demean_id = 'nodemean';
n_regions = 6;

memory_types = {'coarse'}; % 'precise', 'coarse'

res_dat = load(fullfile(res_dir, ...
    ['res_mantel_mht-word_rps-spc-word_subjs206_' demean_id '.mat']));

[~, out_name, ~] = fileparts(fullfile(res_dir, ...
    ['res_mantel_mht-word_rps-spc-word_subjs206_' demean_id '.mat']));

for imem = 1:length(memory_types)
    memory_id = memory_types{imem};
    
    res_r = atanh(res_dat.(['pearsonr_' memory_id]));
    res_p = res_dat.(['p_' memory_id]);
    
    % get the original cluster
    corr_Hmap = res_p <= 0.05;
    
    % extract significant clusters
    for iregion = 1:n_regions
        [L, ~] = bwlabeln(corr_Hmap(iregion,:), 8);
        cluster_value = unique(L);
        cluster_value = setdiff(cluster_value, 0);
        if isempty(cluster_value)   % No significant clusters
            Corr_BehNeu_Orig.(memory_id)(iregion).clusters_idx = [];
            Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu = [];
        else
            % calculate the original r-sum for significant clusters
            for ic = 1:length(cluster_value)
                [i] = find(L == cluster_value(ic));
                r_orig_sum_tmp = res_r(iregion,i);
                % get x-axis & y-axis of the significant areas
                Corr_BehNeu_Orig.(memory_id)(iregion).clusters_idx{ic,1} = i;
                Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu(ic) = sum(r_orig_sum_tmp);
                clear i r_orig_sum_tmp
            end
        end
        clear ic L num cluster_value
    end
    
    % Compare original_sig_clusters with the perm_max_cluster to calculate p-values
    
    % load r sum max
    
    perm_files = dir(fullfile(res_dir, ['mantel_perm_' ...
            memory_id '_multiprcs'], 'res_mantel_mht-word-coarse_isps-word_subjs206_nodemean_z*.mat'));
        n_filess = length(perm_files);
        if n_filess == 1
            load(fullfile(perm_files(1).folder, perm_files(1).name), ...
                'rsum_max_mht')
            rmax_perm_all = rsum_max_mht;
        else
            rmax_perm_all = zeros(n_filess, 1);
            for iperm = 1:n_filess
                load(fullfile(perm_files(iperm).folder, perm_files(iperm).name), ...
                    'rsum_max_mht')
                rmax_perm_all(iperm) = rsum_max_mht;
            end
        end
        n_perms = length(rmax_perm_all);
    
    % R cluster
    for iregion = 1:n_regions
        if ~isempty(Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu)
            for ic = 1:length(Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu)
                Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster(ic) = length(find(rmax_perm_all > Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu(ic)))./n_perms;
            end
        else
            Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster = [];
            Corr_BehNeu_Orig.(memory_id)(iregion).p_min_cluster = [];
        end
        
        Corr_BehNeu_Orig.(memory_id)(iregion).p_max_idx = find(Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster <= 0.05);
        Corr_BehNeu_Orig.(memory_id)(iregion).P_cluster_idx = [Corr_BehNeu_Orig.(memory_id)(iregion).p_max_idx];
    end
    
end
save(fullfile(res_dir, ['zcluster_' out_name '.mat']), 'Corr_BehNeu_Orig')


%% vis only coarse memory region 5 marginal significant (random selected 10 times)
clear,clc;
bids_dir = 'D:\Workplace\shared_memory';
res_dir = fullfile(bids_dir, 'code', '05_stats_semantic_memory', 'results');

addpath(fullfile(bids_dir, 'code', '05_stats_semantic_memory')) % add circle_cluster.m function

demean_id = 'nodemean';
n_regions = 6;

memory_types = {'precise', 'coarse'}; % 'precise', 'coarse'

load(fullfile(res_dir, ...
    ['res_mantel_mht-word_rps-spc-word_perm10mean_subjs206_' demean_id '.mat']), ...
    'pearson_r', 'pvalue');

load(fullfile(res_dir, ...
    ['cluster_res_mantel_mht-word_rps-spc-word_perm10mean_subjs206_' demean_id '.mat']), ...
    'Corr_BehNeu_Orig');

figure('Color','white');
set(gcf,'position',[0 0 680 640]); % set figure position & scale

for iconf = 1:length(memory_types)
    
    subplot(2,1,iconf); % plot subfigure
    
    memory_id = memory_types{iconf};
    
    res_r = pearson_r.(memory_id);
    res_p = pvalue.(memory_id);
    
    i_all = []; j_all = [];
    for iregion = 1:n_regions
        p_tmp = Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster;
        if ~isempty(p_tmp)
            for ip = 1:length(p_tmp)
                if p_tmp(ip) <= 0.1
                    j = Corr_BehNeu_Orig.(memory_id)(iregion).clusters_idx{ip, 1};
                    j_all = [j_all j];
                    i = repmat(iregion, 1, length(j));
                    i_all = [i_all i];
                end
            end
        end
    end
    imagesc(res_r);axis xy;colormap('jet');%caxis([-0.03 0.03]);
    colorbar
    xticks([1 11 21.3 31.5 41.8]);
    set(gca,'xticklabel',[0 200 400 600 800]);
    yticks(1:6);   % y axis position
    set(gca,'yticklabel',1:6);
    xlabel('Encoding time (ms)'), ylabel('Region')
    set(gca,'linewidth',1.5,'fontsize',16,'fontname','arial','ticklength',[0.01 0.02]);
    title([memory_id ' cluster-sig']);
    hold on;
    circle_cluster(j_all,i_all,'k',1.5);
end

%% vis only coarse memory region 5 marginal significant (206 raw data)
clear,clc;
bids_dir = 'D:\Workplace\shared_memory';
res_dir = fullfile(bids_dir, 'code', '05_stats_semantic_memory', 'results');

addpath(fullfile(bids_dir, 'code', '05_stats_semantic_memory')) % add circle_cluster.m function

demean_id = 'nodemean';
n_regions = 6;

memory_types = {'precise', 'coarse'}; % 'precise', 'coarse'

res_dat = load(fullfile(res_dir, ...
    ['res_mantel_mht-word_rps-spc-word_subjs206_' demean_id '.mat']));

load(fullfile(res_dir, ...
    ['cluster_res_mantel_mht-word_rps-spc-word_subjs206_' demean_id '.mat']), ...
    'Corr_BehNeu_Orig');

figure('Color','white');
set(gcf,'position',[0 0 700 720]); % set figure position & scale

for iconf = 1:length(memory_types)
    
    subplot(2,1,iconf); % plot subfigure
    
    memory_id = memory_types{iconf};
    
    res_r = res_dat.(['pearsonr_' memory_id]);
    res_p = res_dat.(['p_' memory_id]);
    
    i_all = []; j_all = [];
    for iregion = 1:n_regions
        p_tmp = Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster;
        if ~isempty(p_tmp)
            for ip = 1:length(p_tmp)
                if p_tmp(ip) <= 0.1
                    j = Corr_BehNeu_Orig.(memory_id)(iregion).clusters_idx{ip, 1};
                    j_all = [j_all j];
                    i = repmat(iregion, 1, length(j));
                    i_all = [i_all i];
                end
            end
        end
    end
    imagesc(res_r);axis xy;colormap('jet');caxis([-0.03 0.03]);
    colorbar
    xticks([1 11 21.3 31.5 41.8]);
    set(gca,'xticklabel',[0 200 400 600 800]);
    yticks(1:6);   % y axis position
    set(gca,'yticklabel',1:6);
    xlabel('Encoding time (ms)'), ylabel('Region')
    set(gca,'linewidth',2,'fontsize',22,'fontname','Gill Sans MT','ticklength',[0.01 0.02]);
    title([memory_id ' cluster-sig']);
    hold on;
    circle_cluster(j_all,i_all,'k',2);
end


%% vis partial mantel memory and representational space (186 raw data)
clear,clc;
bids_dir = 'D:\Workplace\shared_memory\code\05_stats_semantic_memory';
res_dir = fullfile(bids_dir, 'results');

addpath(fullfile(bids_dir)) % add circle_cluster.m function

demean_id = 'nodemean';
n_regions = 6;

memory_types = {'coarse'}; % 'precise', 'coarse'
conf_types = {'dprime_rk'}; % 'age', 'sex', 

for imem = 1:length(memory_types)
    memory_id = memory_types{imem};
    
    load(fullfile(res_dir, ...
        ['res_parmantel_dprime_mht-word-' memory_id '_isps-word_subjs206_' ...
        demean_id '1w.mat']));
    
    [~, out_name, ~] = fileparts(fullfile(res_dir, ...
        ['res_parmantel_dprime_mht-word-' memory_id '_isps-word_subjs206_' ...
        demean_id '1w.mat']));
    
    for iconf = 1:length(conf_types)
        conf_id = conf_types{iconf};
        res_r = atanh(pearson_r.([conf_id '_dist'])); %%% Fisher's Z
        res_p = pvalue.([conf_id '_dist']);
        
        % get the original cluster
        corr_Hmap = res_p <= 0.05;
        
        % extract significant clusters
        for iregion = 1:n_regions
            [L, ~] = bwlabeln(corr_Hmap(iregion,:), 8);
            cluster_value = unique(L);
            cluster_value = setdiff(cluster_value, 0);
            if isempty(cluster_value)   % No significant clusters
                Corr_BehNeu_Orig.(conf_id)(iregion).clusters_idx = [];
                Corr_BehNeu_Orig.(conf_id)(iregion).r_orig_sum_sigclu = [];
            else
                % calculate the original r-sum for significant clusters
                for ic = 1:length(cluster_value)
                    [i] = find(L == cluster_value(ic));
                    r_orig_sum_tmp = res_r(iregion,i);
                    % get x-axis & y-axis of the significant areas
                    Corr_BehNeu_Orig.(conf_id)(iregion).clusters_idx{ic,1} = i;
                    Corr_BehNeu_Orig.(conf_id)(iregion).r_orig_sum_sigclu(ic) = sum(r_orig_sum_tmp);
                    clear i r_orig_sum_tmp
                end
            end
            clear ic L num cluster_value
        end
        
        
        % Compare original_sig_clusters with the perm_max_cluster to calculate p-values
        
        
        % load r sum max
        
        perm_files = dir(fullfile(res_dir, ['parmantel_dprime_perm_' ...
            memory_id '_multiprcs'], 'res_parmantel_dprime_mht-word-coarse_isps-word_subjs206_nodemean_z*.mat'));
        n_filess = length(perm_files);
        if n_filess == 1
            load(fullfile(perm_files(1).folder, perm_files(1).name), ...
                'rsum_max_mht')
            rmax_perm_all = rsum_max_mht.([conf_id '_dist'])';
        else
            rmax_perm_all = zeros(n_filess, 1);
            for iperm = 1:n_filess
                load(fullfile(perm_files(iperm).folder, perm_files(iperm).name), ...
                    'rsum_max_mht')
                rmax_perm_all(iperm) = rsum_max_mht.([conf_id '_dist']);
            end
        end
        n_perms = length(rmax_perm_all);
        
        % R cluster
        for iregion = 1:n_regions
            if ~isempty(Corr_BehNeu_Orig.(conf_id)(iregion).r_orig_sum_sigclu)
                for ic = 1:length(Corr_BehNeu_Orig.(conf_id)(iregion).r_orig_sum_sigclu)
                    Corr_BehNeu_Orig.(conf_id)(iregion).p_max_cluster(ic) = length(find(rmax_perm_all > Corr_BehNeu_Orig.(conf_id)(iregion).r_orig_sum_sigclu(ic)))./n_perms;
                end
            else
                Corr_BehNeu_Orig.(conf_id)(iregion).p_max_cluster = [];
                Corr_BehNeu_Orig.(conf_id)(iregion).p_min_cluster = [];
            end
            
            Corr_BehNeu_Orig.(conf_id)(iregion).p_max_idx = find(Corr_BehNeu_Orig.(conf_id)(iregion).p_max_cluster <= 0.05);
            Corr_BehNeu_Orig.(conf_id)(iregion).P_cluster_idx = [Corr_BehNeu_Orig.(conf_id)(iregion).p_max_idx];
        end
        
    end
    save(fullfile(res_dir, ['zcluster_' out_name '.mat']), 'Corr_BehNeu_Orig')
end


%% vis only coarse memory region 5 marginal significant (partial mantel)
%clear,clc;
bids_dir = 'D:\Workplace\shared_memory\code\05_stats_semantic_memory';
res_dir = fullfile(bids_dir, 'results');

addpath(fullfile(bids_dir)) % add circle_cluster.m function

demean_id = 'nodemean';
n_regions = 6;

memory_id = 'coarse'; % 'precise', 'coarse'

load(fullfile(res_dir, ...
    ['res_parmantel_dprime_mht-word-' memory_id '_isps-word_subjs206_' ...
    demean_id '1w.mat']));

conf_types = {'dprime_rk'};

load(fullfile(res_dir, ...
    ['zcluster_res_parmantel_dprime_mht-word-' memory_id '_isps-word_subjs206_' ...
    demean_id '1w.mat']), 'Corr_BehNeu_Orig');

figure('Color','white');
set(gcf,'position',[0 0 680 640]); % set figure position & scale

for iconf = 1:length(conf_types)
    
    subplot(2,1,iconf); % plot subfigure
    
    conf_id = conf_types{iconf};
    res_r = pearson_r.([conf_id '_dist']);
    res_p = pvalue.([conf_id '_dist']);
    
    i_all = []; j_all = [];
    for iregion = 1:n_regions
        p_tmp = Corr_BehNeu_Orig.(conf_id)(iregion).p_max_cluster;
        if ~isempty(p_tmp)
            for ip = 1:length(p_tmp)
                if p_tmp(ip) <= 0.1
                    j = Corr_BehNeu_Orig.(conf_id)(iregion).clusters_idx{ip, 1};
                    j_all = [j_all j];
                    i = repmat(iregion, 1, length(j));
                    i_all = [i_all i];
                end
            end
        end
    end
    imagesc(res_r);axis xy;colormap('jet');%caxis([-0.03 0.03]);
    colorbar
    xticks([1 11 21.3 31.5 41.8]);
    set(gca,'xticklabel',[0 200 400 600 800]);
    yticks(1:6);   % y axis position
    set(gca,'yticklabel',1:6);
    xlabel('Encoding time (ms)'), ylabel('Region')
    set(gca,'linewidth',1.5,'fontsize',16,'fontname','arial','ticklength',[0.01 0.02]);
    title([conf_id ' cluster-sig']);
    hold on;
    circle_cluster(j_all,i_all,'k',1.5);
end

%% vis hitrate/dprime and individual_group_mean_corr (206subjs)

clear;
bids_dir = 'D:\Workplace\shared_memory';
res_dir = fullfile(bids_dir, 'code', '05_stats_semantic_memory', 'results');

addpath(fullfile(bids_dir, 'code', '05_stats_semantic_memory')) % add circle_cluster.m function

demean_id = 'nodemean';

memory_types = {'coarse'}; % 'precise', 'coarse'

res_dat = load(fullfile(res_dir, ...
    ['corr_grp_sub_rps_spc_dprime_subjs206_' demean_id '.mat']));

n_regions = size(res_dat.R_real, 1);

memory_id = memory_types{1};

res_r = res_dat.R_real;
res_p = res_dat.P_real;

% get the original cluster
corr_Hmap = res_p <= 0.05;

% extract significant clusters
for iregion = 1:n_regions
    [L, ~] = bwlabeln(corr_Hmap(iregion,:), 8);
    cluster_value = unique(L);
    cluster_value = setdiff(cluster_value, 0);
    if isempty(cluster_value)   % No significant clusters
        Corr_BehNeu_Orig.(memory_id)(iregion).clusters_idx = [];
        Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu = [];
    else
        % calculate the original r-sum for significant clusters
        for ic = 1:length(cluster_value)
            [i] = find(L == cluster_value(ic));
            r_orig_sum_tmp = res_r(iregion,i);
            % get x-axis & y-axis of the significant areas
            Corr_BehNeu_Orig.(memory_id)(iregion).clusters_idx{ic,1} = i;
            Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu(ic) = sum(r_orig_sum_tmp);
            clear i r_orig_sum_tmp
        end
    end
    clear ic L num cluster_value
end

% Compare original_sig_clusters with the perm_max_cluster to calculate p-values

% load r sum max
rmax_perm_all = res_dat.rsum_max_mean;
n_perms = 1000;
% R cluster
for iregion = 1:n_regions
    if ~isempty(Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu)
        for ic = 1:length(Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu)
            Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster(ic) = length(find(rmax_perm_all > Corr_BehNeu_Orig.(memory_id)(iregion).r_orig_sum_sigclu(ic)))./n_perms;
        end
    else
        Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster = [];
        Corr_BehNeu_Orig.(memory_id)(iregion).p_min_cluster = [];
    end
    
    Corr_BehNeu_Orig.(memory_id)(iregion).p_max_idx = find(Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster <= 0.05);
    Corr_BehNeu_Orig.(memory_id)(iregion).P_cluster_idx = [Corr_BehNeu_Orig.(memory_id)(iregion).p_max_idx];
end

figure('Color','white');
set(gcf,'position',[0 0 680 640]); % set figure position & scale

for iconf = 1:1
    
    subplot(2,1,iconf); % plot subfigure
    
    i_all = []; j_all = [];
    for iregion = 1:n_regions
        p_tmp = Corr_BehNeu_Orig.(memory_id)(iregion).p_max_cluster;
        if ~isempty(p_tmp)
            for ip = 1:length(p_tmp)
                if p_tmp(ip) <= 0.05
                    j = Corr_BehNeu_Orig.(memory_id)(iregion).clusters_idx{ip, 1};
                    j_all = [j_all j];
                    i = repmat(iregion, 1, length(j));
                    i_all = [i_all i];
                end
            end
        end
    end
    imagesc(res_r);axis xy;colormap('jet');%caxis([-0.03 0.03]);
    colorbar
    xticks([1 11.8 22 32.4 42.8]);
    set(gca,'xticklabel',[0 200 400 600 800]);
    yticks(1:6);   % y axis position
    set(gca,'yticklabel',1:6);
    xlabel('Encoding time (ms)'), ylabel('Region')
    set(gca,'linewidth',2,'fontsize',20,'fontname','Gill Sans MT','ticklength',[0.01 0.02]);
    %title([conf_id ' cluster-sig']);
    hold on;
    circle_cluster(j_all,i_all,'k',2);
end
