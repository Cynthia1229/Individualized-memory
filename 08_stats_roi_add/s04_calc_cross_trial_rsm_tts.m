%%%% merge matrix to group
roi_names = {'PMC', 'lHIP', 'rHIP', 'vmPFC'};

bids_dir = '/seastor/jintaosheng1229/shared_memory';
deriv_dir = fullfile(bids_dir, 'derivatives');

data_dir = fullfile(deriv_dir, 'within-item_cross-subject_PS', 'group');
grp_behav = readtable(fullfile(data_dir, 'group_behav_labels_subj415.csv'));

sig_dir = fullfile(deriv_dir, 'parameters_singletrial-leastsquare-separate', ...
    'signal_roi-add_nodemean', 'group');

subj_all = unique(grp_behav.subj_id);
n_subjs = length(subj_all);
[subj_id1, subj_id2] = meshgrid(subj_all);

tril_idx = logical(tril(ones(n_subjs), -1));
n_pairs = sum(tril_idx(:));

face_list = unique(grp_behav.face_id);
n_faces = length(face_list);

out_tbl = [];
for i_roi = 1:length(roi_names)
    roi_name = roi_names{i_roi};
    fprintf('Merging data for %s ... \n', roi_name)
    load(fullfile(sig_dir, ['group_ste_signal_' roi_name '_subjs415.mat']), 'grp_sig_mat')
    grp_tts = squeeze(nanmean(grp_sig_mat, 2)); %#ok<NANMEAN>
    n_trials = size(grp_sig_mat, 1);
    tril_idx_ct = logical(tril(ones(n_trials), -1));
    grp_ct_tril = zeros(sum(tril_idx_ct(:)), n_subjs);
    for i_sub = 1:n_subjs
        sig_tmp = grp_sig_mat(:,:,i_sub);
        ct_rsm = corr(sig_tmp', 'rows', 'complete');
        grp_ct_tril(:,i_sub) = ct_rsm(tril_idx_ct);
    end
    save(fullfile(bids_dir, 'code', '08_stats_roi_add', 'data', ...
        ['grp_rsm_tts_' roi_name '_subjs415.mat']), 'grp_tts', 'grp_ct_tril')

    isr_ct_rsm = corr(atanh(grp_ct_tril), 'rows', 'complete');
    isr_tts = corr(grp_tts, 'rows', 'complete');

    tmp_tbl = table();
    tmp_tbl.subj_id1 = subj_id1(tril_idx);
    tmp_tbl.subj_id2 = subj_id2(tril_idx);
    tmp_tbl.roi_name = repmat({roi_name}, n_pairs, 1);
    tmp_tbl.gisps = isr_ct_rsm(tril_idx);
    tmp_tbl.isc = isr_tts(tril_idx);
    out_tbl = [out_tbl; tmp_tbl]; %#ok<AGROW>
end
writetable(out_tbl, fullfile(bids_dir, 'code', '08_stats_roi_add', 'data', ...
    'isr_gisps_isc_subjs415.csv'))



