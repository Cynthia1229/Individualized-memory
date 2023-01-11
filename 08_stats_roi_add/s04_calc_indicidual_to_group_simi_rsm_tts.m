%%%% merge matrix to group
roi_names = {'PMC', 'lHIP', 'rHIP', 'vmPFC'};

bids_dir = '/seastor/jintaosheng1229/shared_memory';
subjs_behav = readtable(fullfile(bids_dir, 'facename_behav.csv'));
subj_list = subjs_behav.subj_id(strcmp(subjs_behav.lag, '8_15'),:);

out_tbl = table;

for i_roi = 1:length(roi_names)
    roi_name = roi_names{i_roi};

    load(fullfile(bids_dir, 'code', '08_stats_roi_add', 'data', ...
        ['grp_rsm_tts_' roi_name '_subjs415.mat']), 'grp_tts', 'grp_ct_tril')

    grp_ct_tril_fz = atanh(grp_ct_tril);
    grsm = mean(grp_ct_tril_fz, 2);
    corr_indiv_grsm = corr(grp_ct_tril_fz, grsm);

    gtts = mean(grp_tts, 2);
    corr_indiv_gtts = corr(grp_tts, gtts);

    tbl_tmp = table;
    tbl_tmp.subj_id = subj_list;
    tbl_tmp.roi_name = repmat({roi_name}, length(subj_list), 1);
    tbl_tmp.corr_indiv_grsm = corr_indiv_grsm;
    tbl_tmp.corr_indiv_gtts = corr_indiv_gtts;

    out_tbl = [out_tbl; tbl_tmp];
end
writetable(out_tbl, fullfile(bids_dir, 'code', '08_stats_roi_add', 'data', ...
    'corr_individual_to_grp_subjs415.csv'))
