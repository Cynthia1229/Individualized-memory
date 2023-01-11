bids_dir = 'D:\Workplace\shared_memory\code\08_stats_roi_add';
dat_dir = fullfile(bids_dir, 'data');
load(fullfile(dat_dir, 'compare_rrvsff_subjpairs_matrix.mat'))

% perm_num = 100;

idx_tril = logical(tril(ones(size(rr.lHIP)), -1));

rois = fieldnames(ff);
out_tbl = table();
for i_roi = 1:length(rois)
    tmp_tbl = table();
    roi = rois{i_roi};
    diff_tmp = rr.(roi) - ff.(roi);
    diff_mean_real.(roi) = mean(diff_tmp(idx_tril));
    [pval, t_orig, crit_t, est_alpha, seed_state] = ...
        mult_comp_perm_t1(diff_tmp(idx_tril),10000,0,0.05,0,0);
    tmp_tbl.roi = {roi};
    tmp_tbl.p_value = pval;
    tmp_tbl.t_orig = t_orig;
    tmp_tbl.crit_t_lower = crit_t(1);
    tmp_tbl.crit_t_upper = crit_t(2);
    tmp_tbl.est_alpha = est_alpha;
    %tmp_tbl.diff_real = diff_mean_real;
    out_tbl = [out_tbl; tmp_tbl];
end
writetable(out_tbl, fullfile(bids_dir, 'results', 'res_stats_rr_vs_ff_perm1w.csv'))
% 
% for i_perm = 1:perm_num
%     
%     shuffle matrix
%     ipe = randperm(size(rr.lANG, 1));
%     
%     for i_roi = 1:length(rois)
%         roi = rois{i_roi};
%         rr_mat = rr.(roi);
%         ff_mat = ff.(roi);
%         ff_temp = ff_mat(ipe, ipe);
%         diff_perm = rr_mat - ff_temp;
%         diff_mean_perm.(roi)(i_perm) = mean(diff_perm(idx_tril));
%     end
% end
 
% out_tbl = table();
% for i_roi = 1:length(rois)
%     tmp_tbl = table();
%     roi = rois{i_roi};
%     tmp_tbl.roi = {roi};
%     diff_real = diff_mean_real.(roi);
%     diff_perm = diff_mean_perm.(roi);
%     tmp_tbl.p_value = length(find(diff_perm > diff_real))/length(diff_perm);
%     tmp_tbl.ci_upper = prctile(diff_perm, 95, "all");
%     tmp_tbl.ci_lower = prctile(diff_perm, 5, "all");
%     tmp_tbl.diff_real = diff_real;
%     out_tbl = [out_tbl; tmp_tbl];
% end




