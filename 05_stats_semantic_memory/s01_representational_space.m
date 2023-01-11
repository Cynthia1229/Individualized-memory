function s01_representational_space(sub_num, demean_id)

%%%% Represenatational spce %%%%
base_dir = 'D:\Workplace\shared_memory\derivatives';
data_dir = fullfile(base_dir, 'semantic_memory');

% demean_id = 'demean';

%%% loading data
% data structure: channel * time * trial
in_dt = load(fullfile(data_dir, ['sub' num2str(sub_num) '_SM_rawEEG_demeanedEEG.mat']), ...
    'EEG_replace_badt'); 
EEG_replace_badt = in_dt.EEG_replace_badt;
n_trials = size(EEG_replace_badt, 3);
n_times = size(EEG_replace_badt, 2); % 100/26*307 ms
n_regions = 6; % classical parcellation

% get lower triangular matrix index
tril_idx = logical(tril(ones(n_trials, n_trials), -1));

%%% parameters of constructing representational space %%%
window_length = 26; % 100 ms
step = 5; % 20 ms
window_sets = 1:step:(n_times-window_length+1);

region_infos = readtable(fullfile(data_dir, 'eeg_channel_labels_64.csv'));

if strcmp(demean_id, 'demean')
    EEG_replace_badt = EEG_replace_badt - nanmean(EEG_replace_badt, 3);
end

rps_spc_region_time_fisherz = zeros(n_regions, length(window_sets), ...
    n_trials*(n_trials-1)/2);

for iregion = 1:n_regions
    % get region index
    region_idx = region_infos{:, ['region' num2str(iregion)]} == iregion;
    for wd = 1:length(window_sets)
        % get time window idx
        time_idx = window_sets(wd):(window_sets(wd)+window_length-1);
        dat_tmp = reshape(EEG_replace_badt(region_idx, time_idx, :), ...
            [sum(region_idx) * length(time_idx), n_trials]);
        pearson_r = corr(dat_tmp, 'rows', 'pairwise');
        fisher_zr = atanh(pearson_r(tril_idx)); % fisher's r-to-z
        rps_spc_region_time_fisherz(iregion, wd, :) = fisher_zr;
    end
end

% save out
out_dir = fullfile(data_dir, ['SM_rps_spc_' demean_id]);
if ~ exist(out_dir, 'dir')
    mkdir(out_dir)
end

save(fullfile(out_dir, ['sub-' num2str(sub_num) '_rps_spc_' demean_id '.mat']), ...
    'rps_spc_region_time_fisherz')
return



