%%%% merge matrix to group
roi_names = {'PMC', 'lHIP', 'rHIP', 'vmPFC'};

bids_dir = '/seastor/jintaosheng1229/shared_memory';
deriv_dir = fullfile(bids_dir, 'derivatives');

sig_dir = fullfile(deriv_dir, 'parameters_singletrial-leastsquare-separate', ...
    'signal_roi-add_nodemean');
out_dir = fullfile(sig_dir, 'group');
if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end

grp_behav = readtable(fullfile(bids_dir, 'facename_behav.csv'));
subj_list = grp_behav.subj_id(strcmp(grp_behav.lag, '8_15'),:);

for i_roi = 1:length(roi_names)
    roi_name = roi_names{i_roi};
    fprintf('Merging data for %s ... \n', roi_name)
%     grp_roi_signals = table;
    load(fullfile(sig_dir, 'sub-2001', ['sub-2001_space-T1w-to-MNI_signal-' ...
            roi_name '.mat']), 'roi_signal');
    grp_sig_mat = zeros(size(roi_signal,1), size(roi_signal,2), length(subj_list));
    for i_sub = 1:length(subj_list)
        %subj_behav = grp_behav(grp_behav.subj_id == subj_list(i_sub), :);
        subj_name = ['sub-' num2str(subj_list(i_sub))];
        
        load(fullfile(sig_dir, subj_name, [subj_name '_space-T1w-to-MNI_signal-' ...
            roi_name '.mat']), 'roi_signal');
        grp_sig_mat(:,:,i_sub) = roi_signal;
%         tmp_tbl = table;
%         tmp_tbl.ste_signal = roi_signal(:);
%         tmp_tbl.trial_id = repmat([1:size(roi_signal,1)]', size(roi_signal,2), 1); %#ok<NBRAK>
%         tmp_tbl.voxel_id = repelem([1:size(roi_signal,2)]', size(roi_signal,1)); %#ok<NBRAK>
%         tmp_tbl.subj_id = repmat(subj_list(i_sub), length(roi_signal(:)));
%         tmp_tbl.roi_name = repmat(roi_name, length(roi_signal(:)));
%         
%         grp_roi_signals = [grp_roi_signals; tmp_tbl]; %#ok<AGROW>
    end
    
%     writetable(grp_roi_signals, fullfile(sig_dir, 'group', ...
%         ['group_ste_signal_' roi_name '_subjs415.csv']))

    save(fullfile(sig_dir, 'group', ['group_ste_signal_' roi_name '_subjs415']), ...
        'grp_sig_mat', '-mat')
end
