%%%% correlation between remove ERP and non-remove ERP

base_dir = 'D:\Workplace\shared_memory\derivatives';
data_dir = fullfile(base_dir, 'semantic_memory');

subj_list = [2002, 2005];
n_subjs = length(subj_list);

load(fullfile(data_dir, 'SM_rps_spc_demean', 'sub-2002_rps_spc_demean.mat'), ...
    'rps_spc_region_time_fisherz');
n_regions = size(rps_spc_region_time_fisherz, 1);
n_windows = size(rps_spc_region_time_fisherz, 2);

corr_rps_spc = zeros(n_regions, n_windows, n_subjs);

for isub = 1:length(subj_list)
    % loading data with removing ERP
    dt_rmerp = load(fullfile(data_dir, 'SM_rps_spc_demean', ...
        ['sub-' num2str(subj_list(isub)) '_rps_spc_demean.mat']), ...
        'rps_spc_region_time_fisherz');
    dt_rmerp = dt_rmerp.rps_spc_region_time_fisherz;
    % loading data without removing ERP
    dt_nrmerp = load(fullfile(data_dir, 'SM_rps_spc_nondemean', ...
        ['sub-' num2str(subj_list(isub)) '_rps_spc_nondemean.mat']), ...
        'rps_spc_region_time_fisherz');
    dt_nrmerp = dt_nrmerp.rps_spc_region_time_fisherz;
    
    for iregion = 1:n_regions
        for iwd = 1:n_windows
            corr_rps_spc(iregion, iwd, isub) = corr(squeeze(dt_rmerp(iregion, iwd, :)), ...
                squeeze(dt_nrmerp(iregion, iwd, :)), 'rows', 'pairwise');
        end
    end
end

