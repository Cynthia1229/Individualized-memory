%%% only calculating Pearson correlation
base_dir = '/seastor/jintaosheng1229/eeg_word_memory';
deriv_dir = fullfile(base_dir, 'derivatives');

subjs_mht = readtable(fullfile(base_dir, 'subj206_pairs_mht_dist.csv'));

memory_types = {'precise', 'coarse'};

dist_rps_spc_demean = load(fullfile(deriv_dir, 'rps_spc_demean', 'group', ...
    'dist_cross_subj_rps_spc_pairwise_demean.mat'), 'dist_cross_subj_rps_spc');

dist_rps_spc_nodemean = load(fullfile(deriv_dir, 'rps_spc_nodemean', 'group', ...
    'dist_cross_subj_rps_spc_pairwise_nodemean.mat'), 'dist_cross_subj_rps_spc');

n_regions = size(dist_rps_spc_demean.dist_cross_subj_rps_spc, 1);
n_times = size(dist_rps_spc_demean.dist_cross_subj_rps_spc, 2);

n_subjs = size(dist_rps_spc_demean.dist_cross_subj_rps_spc, 3);
tril_idx = logical(tril(ones(n_subjs), -1));
n_pairs = sum(tril_idx(:));

dist_demean = zeros(n_regions, n_times, n_pairs);
dist_nodemean = zeros(n_regions, n_times, n_pairs);
for iregion = 1:n_regions
    for itime = 1:n_times
        tmp = squeeze(dist_rps_spc_demean.dist_cross_subj_rps_spc(iregion, itime, :, :));
        dist_demean(iregion, itime, :) = tmp(tril_idx);
        tmp = squeeze(dist_rps_spc_nodemean.dist_cross_subj_rps_spc(iregion, itime, :, :));
        dist_nodemean(iregion, itime, :) = tmp(tril_idx);
    end
end

for imemory = 1:2
    memory_id = memory_types{imemory};
    memory_mht = table2array(subjs_mht(strcmp(subjs_mht.memory_type, memory_id), 'distance'));
    for iregion = 1:6
        corr_rps_spc_mht_demean(iregion, :, imemory) = corr(memory_mht, squeeze(dist_demean(iregion, :, :))');
        corr_rps_spc_mht_nodeman(iregion, :, imemory) = corr(memory_mht, squeeze(dist_nodemean(iregion, :, :))');
    end
end

%% plot out

figure('NumberTitle','off', 'Color','white');
set(gcf,'position',[0 0 600 360]); % set figure position & scale
imagesc(corr_rps_spc_mht_nodeman(:,:,1)); axis xy; colormap('autumn'); %caxis([-0.03 0.03]); % summer autumn
colorbar('FontSize', 16, 'LineWidth', 2)

%xticks([1 12 23.5 35 46 57]);   % x axis position [0-57]
xticks([1 11 21.3 31.5 41.8]);
set(gca,'xticklabel',[0 200 400 600 800]);
%yticks([1:1:6]);   % y axis position
%set(gca,'yticklabel',[1:1:6]);
xlabel('Time (ms)'), ylabel('Region')
set(gca,'linewidth',2,'fontsize',16,'fontname','Gill Sans MT','ticklength',[0.01 0.02]);
title('corr neural-withERP and MHT (precise)');
