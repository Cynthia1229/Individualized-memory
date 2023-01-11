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
    load(fullfile(sig_dir, ['group_ste_signal_' roi_name '_subjs415.mat']), 'grp_sig_mat')
    n_voxels = size(grp_sig_mat, 2);
    
    dat_tmp = [];
    for isub = 1:n_subjs
        dat_tmp = [dat_tmp; grp_sig_mat(:,:,isub)]; %#ok<AGROW>
    end
    
    out_tbl = [];
    for iface = 1:n_faces
        face_id = face_list(iface);
        dat_reps = zeros(n_subjs, n_voxels, 2);
        for irep = 1:2
            idx_face = grp_behav.face_id == face_id & grp_behav.num_pres == irep;
            dat_reps(:,:,irep) = dat_tmp(idx_face,:);
        end
        dat_reps_mean = nanmean(dat_reps, 3); %#ok<NANMEAN> 
        corr_wics = corr(dat_reps_mean', 'rows', 'complete');
        tmp_tbl = table();
        tmp_tbl.subj_id1 = subj_id1(tril_idx);
        tmp_tbl.subj_id2 = subj_id2(tril_idx);
        tmp_tbl.face_id = repmat(face_id, n_pairs, 1);
        tmp_tbl.repetition = repmat('mean', n_pairs, 1);
        tmp_tbl.similarity = corr_wics(tril_idx);
        out_tbl = [out_tbl; tmp_tbl]; %#ok<AGROW> 
    end
    writetable(out_tbl, fullfile(data_dir, ...
        ['subjpairs_correlation_within-item_cross-subj_allitems_' roi_name '.csv']))
end
