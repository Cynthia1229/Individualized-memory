%% stats and plot for RSA (based on ERP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/seastor/leoneshi/MEG/circle_cluster.m');
%% predefined parameters
sourcefolder='/seastor/leoneshi/MEG/MEG_analysis/RSA/ERP';
outputfolder='/seastor/leoneshi/MEG/MEG_analysis/RSA/ERP/group/Item_specific';
mkdir(outputfolder);
subject_list=[1:29];
task_list={'semantic','modality'};

%% Visual predictive item-specific
for it=2:length(task_list)
    if it==2
        subject_list=[1:6,8:29]
    end
    for isub=1:length(subject_list)
        subID=subject_list(isub)
        % load WI_Visual data
        load(fullfile(sourcefolder,['sub' num2str(subID,'%02d') '_' task_list{it} '_ERP_WI_Visual_similarity.mat']),'Z_Simi_WI_Visual');
        % avergae across trials within subject
%         for i=1:length(Z_Simi_WI_Visual)
%             simi_tmp1(i,:,:)=Z_Simi_WI_Visual{i};
%         end
        WI_Visual_avg(subID,:,:)=squeeze(mean(Z_Simi_WI_Visual,1));
        clear Z_Simi_WI_Visual

        % load BI_Visual data
        load(fullfile(sourcefolder,['sub' num2str(subID,'%02d') '_' task_list{it} '_ERP_BI_Visual_similarity.mat']),'Z_Simi_BI_Visual');
        % avergae across trials within subject
%         for i=1:length(Z_Simi_BI_Visual)
%             simi_tmp2(i,:,:)=Z_Simi_BI_Visual{i};
%         end
        BI_Visual_avg(subID,:,:)=squeeze(mean(Z_Simi_BI_Visual,1));
        clear Z_Simi_BI_Visual
    end
    save(fullfile(outputfolder,[task_list{it} '_ERP_WI_Visual_similarity_group.mat']),'WI_Visual_avg');
    save(fullfile(outputfolder,[task_list{it} '_ERP_BI_Visual_similarity_group.mat']),'BI_Visual_avg');
    
    %% imagesc plot for each condition and stats
    % imagesc plot for raw similarity
    figure;
    set(gcf,'color','w');
    set(gcf,'position',[0 0 1200 250]);
    subplot(1,3,1);
    imagesc(squeeze(mean(WI_Visual_avg,1)));axis xy;colormap('jet');
%     caxis([0 0.07]);
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'xticklabel',[0.5 1 1.5 2]);
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'yticklabel',[0.5 1 1.5 2]);
    xlabel('1st Run: Item 1 (s)');
    ylabel('2nd Run: Item 1 (s)');
    set(gca,'linewidth',1.5,'fontsize',10,'fontname','Arial','ticklength',[0.01 0.02]);
    title('WI Similarity');
    colorbar;
    
    subplot(1,3,2);
    imagesc(squeeze(mean(BI_Visual_avg,1)));axis xy;colormap('jet');
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'xticklabel',[0.5 1 1.5 2]);
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'yticklabel',[0.5 1 1.5 2]);
    xlabel('1st Run: Item 1 (s)');
    ylabel('2nd Run: Item 2 (s)');
    set(gca,'linewidth',1.5,'fontsize',10,'fontname','Arial','ticklength',[0.01 0.02]);
    title('BI Similarity');
    colorbar;

    % statistical analysis
    [h_real,p_real,~,stats_real]=ttest(WI_Visual_avg,BI_Visual_avg);
    h_map=squeeze(h_real);
    stats_map=squeeze(stats_real.tstat);
    
    subplot(1,3,3);
    imagesc(stats_map);axis xy;colormap('jet');
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'xticklabel',[0.5 1 1.5 2]);
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'yticklabel',[0.5 1 1.5 2]);
    xlabel('Item (s)');
    ylabel('Item (s)');
    set(gca,'linewidth',1.5,'fontsize',10,'fontname','arial','ticklength',[0.01 0.02]);
    title('WI-BI Similarity');
    cb = colorbar;
    hylabel = ylabel(cb,{'T-value'});
    
    % find sig cluster
    [L,num]=bwlabeln(h_map,8);
    cluster_value = unique(L);
    cluster_value = setdiff(cluster_value,0);
    if length(cluster_value)>0
        for ic = 1:length(cluster_value)
            [i,j] = find(L == cluster_value(ic));
            for il = 1:length(i)
                t_orig_sum_tmp(il)= stats_map(i(il),j(il));
            end

            clusters_idx{ic,1} = [i,j];
            clear i j
            t_orig_sum_cluster(ic) = sum(t_orig_sum_tmp);
            clear t_orig_sum_tmp
        end
    end
    
    % permute test
    shuffle_num = 1000;
    data_tmp(:,:,:,1)=WI_Visual_avg;
    data_tmp(:,:,:,2)=BI_Visual_avg;
    clear WI_Visual_avg BI_Visual_avg
    for ishuffle = 1:shuffle_num
        for isub = 1:size(data_tmp,1)
            data_shuffle(isub,:,:,:) = squeeze(data_tmp(isub,:,:,randperm(2)));
        end
        [h,p,~,stats] = ttest(squeeze(data_shuffle(:,:,:,1)),squeeze(data_shuffle(:,:,:,2)));

        H_shuffle = squeeze(h);
        stats_shuffle = squeeze(stats.tstat);
        [L,num]=bwlabeln(H_shuffle,8);
        cluster_value = unique(L);
        cluster_value = setdiff(cluster_value,0);
        if length(cluster_value)>0
            for ic = 1:length(cluster_value)
                [i,j] = find(L == cluster_value(ic));
                cluster_size_tmp = length(i);
                cluster_num(ic) = cluster_size_tmp;

                stats_ttmp = 0;
                for is = 1:cluster_size_tmp
                    stats_tmp = stats_shuffle(i(is),j(is));
                    stats_ttmp = stats_ttmp+stats_tmp;
                    clear stats_tmp
                end
                cluster_stats(ic) = stats_ttmp;
                clear stats_ttmp cluster_size_tmp
            end

            tsum_max_dis(ishuffle) = max(cluster_stats);
            tsum_min_dis(ishuffle) = min(cluster_stats);
            clear cluster_stats cluster_value L num 
        else
            tsum_max_dis(ishuffle) = 0;
            tsum_min_dis(ishuffle) = 0;
        end

        clear sub_idx h p stats H_shuffle stats_shuffle L num cluster_value cluster_stats data_shuffle
    end
    
    for ic = 1:length(t_orig_sum_cluster)
        p_max_cluster(ic)  = length(find(tsum_max_dis>t_orig_sum_cluster(ic)))./shuffle_num;
        p_min_cluster(ic)  = length(find(tsum_min_dis<t_orig_sum_cluster(ic)))./shuffle_num;
    end
    
    p_max_idx = find(p_max_cluster < 0.05);
    p_min_idx = find(p_min_cluster < 0.05);
    P_cluster_idx = [p_max_idx p_min_idx];
    
    hold on;
    for il = 1:length(P_cluster_idx)
        circle_cluster(clusters_idx{P_cluster_idx(il),1}(:,2),clusters_idx{P_cluster_idx(il),1}(:,1),'k',1.5);
        hold on; 
    end
    saveas(gcf,fullfile(outputfolder,[task_list{it} '_Item_specific_Visual.jpg']));
    save(fullfile(outputfolder,[task_list{it} '_Item_specific_Visual_stat.mat']),'h_real','p_real','stats_real','P_cluster_idx','clusters_idx','t_orig_sum_cluster');
    clear p_max_cluster p_min_cluster P_cluster_idx p_max_idx p_min_idx clusters_idx t_orig_sum_cluster data_tmp h_real p_real stats_real
end



%% Auditory item-specific
for it=1:length(task_list)
    if it==2
        subject_list=[1:6,8:29]
    end
    for isub=7:length(subject_list)
        subID=subject_list(isub)
        % load WI_Auditory data
        load(fullfile(sourcefolder,['sub' num2str(subID,'%02d') '_' task_list{it} '_ERP_WI_Auditory_similarity.mat']),'Z_Simi_WI_Auditory');
        % avergae across trials within subject
%         for i=1:length(Z_Simi_WI_Auditory)
%             simi_tmp3(i,:,:)=Z_Simi_WI_Auditory{i};
%         end
        WI_Auditory_avg(subID,:,:)=squeeze(mean(Z_Simi_WI_Auditory,1));
        clear Z_Simi_WI_Auditory

        % load BI_Visual data
        load(fullfile(sourcefolder,['sub' num2str(subID,'%02d') '_' task_list{it} '_ERP_BI_Auditory_similarity.mat']),'Z_Simi_BI_Auditory');
        % avergae across trials within subject
%         for i=1:length(Z_Simi_BI_Auditory)
%             simi_tmp4(i,:,:)=Z_Simi_BI_Auditory{i};
%         end
        BI_Auditory_avg(subID,:,:)=squeeze(mean(Z_Simi_BI_Auditory,1));
        clear Z_Simi_BI_Auditory
    end
    save(fullfile(outputfolder,[task_list{it} '_ERP_WI_Auditory_similarity_group.mat']),'WI_Auditory_avg');
    save(fullfile(outputfolder,[task_list{it} '_ERP_BI_Auditory_similarity_group.mat']),'BI_Auditory_avg');
    
    %% imagesc plot for each condition and stats
    % imagesc plot for raw similarity
    figure;
    set(gcf,'color','w');
    set(gcf,'position',[0 0 1200 250]);
    subplot(1,3,1);
    imagesc(squeeze(mean(WI_Auditory_avg,1)));axis xy;colormap('jet');
%     caxis([0 0.07]);
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'xticklabel',[0.5 1 1.5 2]);
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'yticklabel',[0.5 1 1.5 2]);
    xlabel('1st Run: Item 1 (s)');
    ylabel('2nd Run: Item 1 (s)');
    set(gca,'linewidth',1.5,'fontsize',10,'fontname','Arial','ticklength',[0.01 0.02]);
    title('WI Similarity');
    colorbar;
    
    subplot(1,3,2);
    imagesc(squeeze(mean(BI_Auditory_avg,1)));axis xy;colormap('jet');
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'xticklabel',[0.5 1 1.5 2]);
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'yticklabel',[0.5 1 1.5 2]);
    xlabel('1st Run: Item 1 (s)');
    ylabel('2nd Run: Item 2 (s)');
    set(gca,'linewidth',1.5,'fontsize',10,'fontname','Arial','ticklength',[0.01 0.02]);
    title('BI Similarity');
    colorbar;

    % statistical analysis
    [h_real,p_real,~,stats_real]=ttest(WI_Auditory_avg,BI_Auditory_avg);
    h_map=squeeze(h_real);
    stats_map=squeeze(stats_real.tstat);
    
    subplot(1,3,3);
    imagesc(stats_map);axis xy;colormap('jet');
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'xticklabel',[0.5 1 1.5 2]);
    set(gca,'ytick',[50 100 150 200]);
    set(gca,'yticklabel',[0.5 1 1.5 2]);
    xlabel('Item (s)');
    ylabel('Item (s)');
    set(gca,'linewidth',1.5,'fontsize',10,'fontname','arial','ticklength',[0.01 0.02]);
    title('WI-BI Similarity');
    cb = colorbar;
    hylabel = ylabel(cb,{'T-value'});
    
    % find sig cluster
    [L,num]=bwlabeln(h_map,8);
    cluster_value = unique(L);
    cluster_value = setdiff(cluster_value,0);
    if length(cluster_value)>0
        for ic = 1:length(cluster_value)
            [i,j] = find(L == cluster_value(ic));
            for il = 1:length(i)
                t_orig_sum_tmp(il)= stats_map(i(il),j(il));
            end

            clusters_idx{ic,1} = [i,j];
            clear i j
            t_orig_sum_cluster(ic) = sum(t_orig_sum_tmp);
            clear t_orig_sum_tmp
        end
    end
    
    % permute test
    shuffle_num = 1000;
    data_tmp(:,:,:,1)=WI_Auditory_avg;
    data_tmp(:,:,:,2)=BI_Auditory_avg;
    clear WI_Auditory_avg BI_Auditory_avg
    for ishuffle = 1:shuffle_num
        for isub = 1:size(data_tmp,1)
            data_shuffle(isub,:,:,:) = squeeze(data_tmp(isub,:,:,randperm(2)));
        end
        [h,p,~,stats] = ttest(squeeze(data_shuffle(:,:,:,1)),squeeze(data_shuffle(:,:,:,2)));

        H_shuffle = squeeze(h);
        stats_shuffle = squeeze(stats.tstat);
        [L,num]=bwlabeln(H_shuffle,8);
        cluster_value = unique(L);
        cluster_value = setdiff(cluster_value,0);
        if length(cluster_value)>0
            for ic = 1:length(cluster_value)
                [i,j] = find(L == cluster_value(ic));
                cluster_size_tmp = length(i);
                cluster_num(ic) = cluster_size_tmp;

                stats_ttmp = 0;
                for is = 1:cluster_size_tmp
                    stats_tmp = stats_shuffle(i(is),j(is));
                    stats_ttmp = stats_ttmp+stats_tmp;
                    clear stats_tmp
                end
                cluster_stats(ic) = stats_ttmp;
                clear stats_ttmp cluster_size_tmp
            end

            tsum_max_dis(ishuffle) = max(cluster_stats);
            tsum_min_dis(ishuffle) = min(cluster_stats);
            clear cluster_stats cluster_value L num 
        else
            tsum_max_dis(ishuffle) = 0;
            tsum_min_dis(ishuffle) = 0;
        end

        clear sub_idx h p stats H_shuffle stats_shuffle L num cluster_value cluster_stats data_shuffle
    end
    
    for ic = 1:length(t_orig_sum_cluster)
        p_max_cluster(ic)  = length(find(tsum_max_dis>t_orig_sum_cluster(ic)))./shuffle_num;
        p_min_cluster(ic)  = length(find(tsum_min_dis<t_orig_sum_cluster(ic)))./shuffle_num;
    end
    
    p_max_idx = find(p_max_cluster < 0.1);
    p_min_idx = find(p_min_cluster < 0.1);
    P_cluster_idx = [p_max_idx p_min_idx];
    
    hold on;
    for il = 1:length(P_cluster_idx)
        circle_cluster(clusters_idx{P_cluster_idx(il),1}(:,2),clusters_idx{P_cluster_idx(il),1}(:,1),'k',1.5);
        hold on;
    end
    saveas(gcf,fullfile(outputfolder,[task_list{it} '_Item_specific_Auditory.jpg']));
    save(fullfile(outputfolder,[task_list{it} '_Item_specific_Auditory_stat.mat']),'h_real','p_real','stats_real','P_cluster_idx','clusters_idx','t_orig_sum_cluster');
    clear p_max_cluster p_min_cluster P_cluster_idx p_max_idx p_min_idx clusters_idx t_orig_sum_cluster data_tmp h_real p_real stats_real
    close all
end
clear all