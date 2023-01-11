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

for i_roi = 1:length(roi_names)
    roi_name = roi_names{i_roi};
    fprintf('Merging data for %s ... \n', roi_name)
    load(fullfile(sig_dir, ['group_ste_signal_' roi_name '_subj415.mat']), 'grp_sig_mat')
    n_trials = size(grp_sig_mat, 1);
    n_voxels = size(grp_sig_mat, 2);
    
    dat_tmp = [];
    for isub = 1:n_subjs
        dat_tmp = [dat_tmp; grp_sig_mat(:,:,isub)]; %#ok<AGROW>
    end
    
    corr_wi_cs_reps = zeros(n_subjs, n_subjs, n_faces, 3); % 1-rep1; 2-rep2; 3-rep_mean
    
    for iface = 1:n_faces
        face_id = face_list(iface);
        dat_reps = zeros(n_subjs, n_voxels, 2);
        for irep = 1:2
            idx_face = grp_behav.face_id == face_id & grp_behav.num_pres == irep;
            dat_for_simi = dat_tmp(idx_face,:);
            corr_wi_cs_reps(:, :, iface, irep) = corr(dat_for_simi', 'rows', 'complete');
            dat_reps(:,:,irep) = dat_for_simi;
        end
        dat_reps_mean = nanmean(dat_reps, 3); %#ok<NANMEAN> 
        corr_wi_cs_reps(:, :, iface, 3) = corr(dat_reps_mean', 'rows', 'complete');
    end
    dist_wi_cs_reps_mean = squeeze(nanmean(corr_wi_cs_reps, 3)); %#ok<NANMEAN> 
    out_tbl = table();
    out_tbl.subj_id1 = subj_id1(tril_idx);
    out_tbl.subj_id2 = subj_id2(tril_idx);
    out_tbl.roi_name = repmat(roi_name, n_pairs, 1);
    dist_tmp = dist_wi_cs_reps_mean(:,:,1);
    out_tbl.rep1 = dist_tmp(tril_idx);
    dist_tmp = dist_wi_cs_reps_mean(:,:,2);
    out_tbl.rep2 = dist_tmp(tril_idx);
    dist_tmp = dist_wi_cs_reps_mean(:,:,3);
    out_tbl.repm = dist_tmp(tril_idx);
    writetable(out_tbl, fullfile(data_dir, ...
        ['group_similarity_within-item_cross-subj_mean-items_' roi_name '.csv']))
end
