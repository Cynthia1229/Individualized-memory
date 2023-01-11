%%%%% RAS analysis based on Normalised TF power
%% add toolbox
addpath('/seastor/leoneshi/fieldtrip-20200922');
ft_defaults;
addpath(genpath('/seastor/leoneshi/GroupICATv4.0b'));
%%
sourcefolder='/seastor/leoneshi/MEG/MEG_analysis/';
load('/seastor/leoneshi/MEG/MEG_analysis/encoding_time_tf.mat');
tf=[2:1:29,30:5:100];

%modality task
% subject_list=[22:29];
% for isub=1:length(subject_list)
%     subID=subject_list(isub)
%     if subID==1
%         run_list=[5,6,7];
%     elseif subID==5
%         run_list=[4,5];
%     elseif mod(subID,2)
%         run_list=[4,5,6];
%     else
%         run_list=[1,2,3];
%     end
%     TF_power=[];
%     all_trials=[];
%     for irun=1:length(run_list)
%         runID=run_list(irun);
%         load([sourcefolder filesep 'TF_preprocess' filesep 'Normalization' filesep 'sub' num2str(subID,'%02d') '_TF_modality_norm_run' num2str(runID) '.mat']);
%         % remove wrong response
%         tf_normalized_pcnt_tmp=tf_normalized_pcnt(find(trial_info(:,7)==1),:,tf,:); %only choose 43 freqs
%         trial_info=trial_info(find(trial_info(:,7)==1),:);
%         trial_info(:,1)=irun; %add test RunNum
%         clear tf_normalized_pcnt
%         %
%         TF_power=cat(1,TF_power,tf_normalized_pcnt_tmp); %trial*channel*frequency*time
%         all_trials=cat(1,all_trials,trial_info);
%         clear tf_normalized_pcnt_tmp
%     end
%     
%         %% redefine trial_info
%     % column_name for trial_info: 1-> test RunNum; 2-> test TrialNum; 3 ->
%     % SeqNum; 4-> ItemPos; 5 -> ConNum_Predictable; 6 -> Stimuli; 7 -> ACC;
%     % 8-> RT
%     all_trials(:,9)=[1:length(all_trials)];  %add 9 -> real trialID
%     %% add real runID
%     load(['/seastor/leoneshi/MEG/' filesep 'Rawdata' filesep 'sub-' num2str(subID,'%02d') filesep 'beh' filesep 'sub-' num2str(subID,'%02d') '_task-modality.mat']);
%     for i=1:length(all_trials)
%         for j=1:length(results)
%             if all_trials(i,2)==results(j,3) & all_trials(i,3)==results(j,4) & all_trials(i,4)==results(j,5) & all_trials(i,5)==results(j,6) & all_trials(i,6)==results(j,7) & all_trials(i,8)==results(j,13)
%                 all_trials(i,10)=results(j,2); % add 10 -> real runID
%             end
%         end
%     end
%     clear results i j
%     
%     %% add other properties
%     load(['/seastor/leoneshi/MEG/condition_mapping.mat']);
%     for i=1:length(all_trials)
%         for j=1:size(separate_condition,1)
%             if all_trials(i,3)==separate_condition.SeqNum(j,1) & all_trials(i,4)==separate_condition.ItemPos(j,1)
%                 all_trials(i,11)=separate_condition.ConNum_Predictive(j,1); % add 11 -> ConNum_Predictive
%                 all_trials(i,12)=separate_condition.Modality(j,1); % add 12 -> Modality
%                 all_trials(i,13)=separate_condition.Predictable(j,1); % add 13 -> Predictable
%                 all_trials(i,14)=separate_condition.Shift(j,1); % add 14 -> Shift
%                 all_trials(i,15)=separate_condition.Predictive(j,1); % add 15 -> Predictive
%                 all_trials(i,16)=separate_condition.Shifting(j,1); % add 16 -> Shifitng
%                 all_trials(i,17)=separate_condition.SeqType(j,1); % add 17 -> SeqType
%             end
%         end
%     end
%     clear separate_condition i j
%     
%     %% create pairwise labels for similarity vectors
%     trial_sum=table();
%     % test RunNum
%     [x,y]=meshgrid(all_trials(:,1),all_trials(:,1));
%     trial_sum.i1_RunNum=icatb_mat2vec(x);
%     trial_sum.i2_RunNum=icatb_mat2vec(y);
%     % real Run ID
%     [x,y]=meshgrid(all_trials(:,10),all_trials(:,10));
%     trial_sum.i1_run_id=icatb_mat2vec(x);
%     trial_sum.i2_run_id=icatb_mat2vec(y);
%     % test TrialNum
%     [x,y]=meshgrid(all_trials(:,2),all_trials(:,2));
%     trial_sum.i1_trialNum=icatb_mat2vec(x);
%     trial_sum.i2_trialNum=icatb_mat2vec(y);
%     % real Trial ID
%     [x,y]=meshgrid(all_trials(:,9),all_trials(:,9));
%     trial_sum.i1_trial_id=icatb_mat2vec(x);
%     trial_sum.i2_trial_id=icatb_mat2vec(y);
%     % RT
%     [x,y]=meshgrid(all_trials(:,8),all_trials(:,8));
%     trial_sum.i1_RT=icatb_mat2vec(x);
%     trial_sum.i2_RT=icatb_mat2vec(y);
%     % SeqType
%     [x,y]=meshgrid(all_trials(:,17),all_trials(:,17));
%     trial_sum.i1_SeqType=icatb_mat2vec(x);
%     trial_sum.i2_SeqType=icatb_mat2vec(y);
%     % SeqNum
%     [x,y]=meshgrid(all_trials(:,3),all_trials(:,3));
%     trial_sum.i1_SeqNum=icatb_mat2vec(x);
%     trial_sum.i2_SeqNum=icatb_mat2vec(y);
%     % ItemPos
%     [x,y]=meshgrid(all_trials(:,4),all_trials(:,4));
%     trial_sum.i1_ItemPos=icatb_mat2vec(x);
%     trial_sum.i2_ItemPos=icatb_mat2vec(y);
%     % Stimuli
%     [x,y]=meshgrid(all_trials(:,6),all_trials(:,6));
%     trial_sum.i1_Stimuli=icatb_mat2vec(x);
%     trial_sum.i2_Stimuli=icatb_mat2vec(y); 
%     % ConNum_Predictavle
%     [x,y]=meshgrid(all_trials(:,5),all_trials(:,5));
%     trial_sum.i1_ConNum_Predictable=icatb_mat2vec(x);
%     trial_sum.i2_ConNum_Predictable=icatb_mat2vec(y);
%     % ConNum_Predictive
%     [x,y]=meshgrid(all_trials(:,11),all_trials(:,11));
%     trial_sum.i1_ConNum_Predictive=icatb_mat2vec(x);
%     trial_sum.i2_ConNum_Predictive=icatb_mat2vec(y);
%     % Modality
%     [x,y]=meshgrid(all_trials(:,12),all_trials(:,12));
%     trial_sum.i1_Modality=icatb_mat2vec(x);
%     trial_sum.i2_Modality=icatb_mat2vec(y);
%     % Predictable
%     [x,y]=meshgrid(all_trials(:,13),all_trials(:,13));
%     trial_sum.i1_Predictable=icatb_mat2vec(x);
%     trial_sum.i2_Predictable=icatb_mat2vec(y);
%     % Shift
%     [x,y]=meshgrid(all_trials(:,14),all_trials(:,14));
%     trial_sum.i1_Shift=icatb_mat2vec(x);
%     trial_sum.i2_Shift=icatb_mat2vec(y);
%     % Predictive
%     [x,y]=meshgrid(all_trials(:,15),all_trials(:,15));
%     trial_sum.i1_Predictive=icatb_mat2vec(x);
%     trial_sum.i2_Predictive=icatb_mat2vec(y);
%     % Shifitng
%     [x,y]=meshgrid(all_trials(:,16),all_trials(:,16));
%     trial_sum.i1_Shifitng=icatb_mat2vec(x);
%     trial_sum.i2_Shifitng=icatb_mat2vec(y);
%     
%     clear all_trials
%     save([sourcefolder filesep 'RSA' filesep 'TF' filesep 'sub' num2str(subID,'%02d') '_modality_RSA_TF.mat'],'TF_power','trial_sum','-v7.3');
%     clear TF_power
% end


%% selection interest of condition
subject_list=[5:6,8:29];
for isub=1:length(subject_list)
    subID=subject_list(isub)
    load([sourcefolder filesep 'RSA' filesep 'TF' filesep 'sub' num2str(subID,'%02d') '_modality_RSA_TF.mat']);
    % item-specific
    trial_sum_tmp=trial_sum(find(trial_sum.i1_run_id~=trial_sum.i2_run_id),:); % cross runs
    clear trial_sum
    trial_sum_tmp=trial_sum_tmp(find((trial_sum_tmp.i1_SeqNum>8) & (trial_sum_tmp.i2_SeqNum>8)),:); % random sequences
    trial_sum_tmp=trial_sum_tmp(find(trial_sum_tmp.i1_Modality==trial_sum_tmp.i2_Modality),:); % same modality
    trial_sum_WI_Visual=trial_sum_tmp(find((trial_sum_tmp.i1_Stimuli==trial_sum_tmp.i2_Stimuli) & (trial_sum_tmp.i1_Modality==1)),:); % same item visual modality
    trial_sum_WI_Auditory=trial_sum_tmp(find((trial_sum_tmp.i1_Stimuli==trial_sum_tmp.i2_Stimuli) & (trial_sum_tmp.i1_Modality==2)),:); % same item auditory modality
    trial_sum_BI_Visual=trial_sum_tmp(find((trial_sum_tmp.i1_Stimuli~=trial_sum_tmp.i2_Stimuli) & (trial_sum_tmp.i1_Modality==1)),:); % diff item visual modality
    trial_sum_BI_Auditory=trial_sum_tmp(find((trial_sum_tmp.i1_Stimuli~=trial_sum_tmp.i2_Stimuli) & (trial_sum_tmp.i1_Modality==2)),:); % diff item auditory modality
    clear trial_sum_tmp
    
    %% RAS based on TF_power
    % sliding-window parameters
     window_length=20;% 200 ms
     step=1; % 10 ms
     window_sets=1:step:size(encoding_time,2)-window_length+1;
    
    % WI_Visual similarity
    for pair=1:size(trial_sum_WI_Visual,1)
        for wd=1:length(window_sets)        
            dat_tmp1(wd,:)=reshape(squeeze(mean(TF_power(trial_sum_WI_Visual.i1_trial_id(pair,1),:,:,[window_sets(wd):window_sets(wd)+window_length-1]),4)),1,306*length(tf));
            dat_tmp2(wd,:)=reshape(squeeze(mean(TF_power(trial_sum_WI_Visual.i2_trial_id(pair,1),:,:,[window_sets(wd):window_sets(wd)+window_length-1]),4)),1,306*length(tf));
        end
        similarity=1-pdist2(dat_tmp1,dat_tmp2,'spearman');
        Z_Simi_WI_Visual(pair,:,:)=0.5*log((1+similarity)./(1-similarity));
        clear data_tmp1 data_tmp2 similarity
    end
    save([sourcefolder filesep 'RSA' filesep 'TF' filesep 'sub' num2str(subID,'%02d') '_modality_TF_WI_Visual_similarity.mat'],'Z_Simi_WI_Visual','trial_sum_WI_Visual','window_sets','-v7.3');
    clear trial_sum_WI_Visual Z_Simi_WI_Visual
    
    % WI_Auditory similarity
    for pair=1:size(trial_sum_WI_Auditory,1)
        for wd=1:length(window_sets)        
            dat_tmp1(wd,:)=reshape(squeeze(mean(TF_power(trial_sum_WI_Auditory.i1_trial_id(pair,1),:,:,[window_sets(wd):window_sets(wd)+window_length-1]),4)),1,306*length(tf));
            dat_tmp2(wd,:)=reshape(squeeze(mean(TF_power(trial_sum_WI_Auditory.i2_trial_id(pair,1),:,:,[window_sets(wd):window_sets(wd)+window_length-1]),4)),1,306*length(tf));
        end
        similarity=1-pdist2(dat_tmp1,dat_tmp2,'spearman');
        Z_Simi_WI_Auditory(pair,:,:)=0.5*log((1+similarity)./(1-similarity));
        clear data_tmp1 data_tmp2 similarity
    end
    save([sourcefolder filesep 'RSA' filesep 'TF' filesep 'sub' num2str(subID,'%02d') '_modality_TF_WI_Auditory_similarity.mat'],'Z_Simi_WI_Auditory','trial_sum_WI_Auditory','window_sets','-v7.3');
    clear trial_sum_WI_Auditory Z_Simi_WI_Auditory 

    % BI_Visual similarity
    for pair=1:size(trial_sum_BI_Visual,1)
        for wd=1:length(window_sets)        
            dat_tmp1(wd,:)=reshape(squeeze(mean(TF_power(trial_sum_BI_Visual.i1_trial_id(pair,1),:,:,[window_sets(wd):window_sets(wd)+window_length-1]),4)),1,306*length(tf));
            dat_tmp2(wd,:)=reshape(squeeze(mean(TF_power(trial_sum_BI_Visual.i2_trial_id(pair,1),:,:,[window_sets(wd):window_sets(wd)+window_length-1]),4)),1,306*length(tf));
        end
        similarity=1-pdist2(dat_tmp1,dat_tmp2,'spearman');
        Z_Simi_BI_Visual(pair,:,:)=0.5*log((1+similarity)./(1-similarity));
        clear data_tmp1 data_tmp2 similarity
    end
    save([sourcefolder filesep 'RSA' filesep 'TF' filesep 'sub' num2str(subID,'%02d') '_modality_TF_BI_Visual_similarity.mat'],'Z_Simi_BI_Visual','trial_sum_BI_Visual','window_sets','-v7.3');
    clear trial_sum_BI_Visual Z_Simi_BI_Visual
    
    % BI_Auditory similarity
    for pair=1:size(trial_sum_BI_Auditory,1)
        for wd=1:length(window_sets)        
            dat_tmp1(wd,:)=reshape(squeeze(mean(TF_power(trial_sum_BI_Auditory.i1_trial_id(pair,1),:,:,[window_sets(wd):window_sets(wd)+window_length-1]),4)),1,306*length(tf));
            dat_tmp2(wd,:)=reshape(squeeze(mean(TF_power(trial_sum_BI_Auditory.i2_trial_id(pair,1),:,:,[window_sets(wd):window_sets(wd)+window_length-1]),4)),1,306*length(tf));
        end
        similarity=1-pdist2(dat_tmp1,dat_tmp2,'spearman');
        Z_Simi_BI_Auditory(pair,:,:)=0.5*log((1+similarity)./(1-similarity));
        clear data_tmp1 data_tmp2 similarity
    end
    save([sourcefolder filesep 'RSA' filesep 'TF' filesep 'sub' num2str(subID,'%02d') '_modality_TF_BI_Auditory_similarity.mat'],'Z_Simi_BI_Auditory','trial_sum_BI_Auditory','window_sets','-v7.3');
    clear trial_sum_BI_Auditory Z_Simi_BI_Auditory
    clear TF_power
end