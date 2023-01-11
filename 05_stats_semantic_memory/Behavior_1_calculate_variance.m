%% calculate behavior performance
% Sisi Wang 09/27/2020 modified from YYL's bev_stat code


%% calculated beh variances
% RT, ACC, d-prime
clear;clc;
input_data_path = 'E:\sisi2020\Experiments\Semantic_memory_EEG\ERP_analysis_sisi\Beha_new\';
output_data_path = 'E:\sisi2020\Experiments\Semantic_memory_EEG\ERP_analysis_sisi\Beha_new\beh_results\';

% load sublist
para_path = 'E:\sisi2020\Experiments\Semantic_memory_EEG\ERP_analysis_sisi\SM_EEG_analysis_codes\';
load([para_path 'sub_all_good.mat']);
sub_list = sub_all_good_231;    clear sub_all_good_231

Result = table;

for s = 1:length(sub_list)
    load([input_data_path 'SM' num2str(sub_list(s)) '_beh.mat']);    
    LSM=SM(4:end-3,:); % learning phase exclude the first&last 3 trials
    clear SM
    load([input_data_path 'TEST' num2str(sub_list(s)) '_beh.mat']);    
    TSM=SM;  % test phase
    clear SM
    
    % define abs/con in TSM
    for i = 1:length(TSM)
        % the second colum is abs/con: abs [1:75 500:543 601:625] 145
        if ismember(TSM(i,2),[1:75])==1 || ismember(TSM(i,2),[500:543])==1 || ismember(TSM(i,2),[601:625])==1
            TSM(i,8)=1;% abs 145
        else
            TSM(i,8)=2;% concrete 155
        end
    end
    
    % TSM colum8 absolute & concrete information
   
    % merge LSM & TSM 
    LSM_tmp = sortrows(LSM,2);
    LSM_tmp = LSM_tmp(:,3:6);
    TSM_tmp = sortrows(TSM,2);
    TSM_tmp(1:150,9:12) = LSM_tmp;% abs/con;resp;rt1;right/wrong
    TSM_tmp(151:300,9:12)=999;% for NaN
    SM = sortrows(TSM_tmp,1);
    
    % delete rt<0.2 trials
    LSM = LSM(LSM(:,5)>0.2,:);
    TSM = TSM(TSM(:,6)>0.2,:);
    SM = SM(SM(:,6)>0.2 & SM(:,11)>0.2,:);
    
    old = SM(SM(:,3)==1,:);
    new = SM(SM(:,3)==2,:);
    Abstract_old = old(old(:,8)==1,:);
    Concrete_old = old(old(:,8)==2,:);
    Abstract_new = new(new(:,8)==1,:);
    Concrete_new = new(new(:,8)==2,:);
    
    Result.subID(s,1) = {sub_list(s)};
    % RT: (Abs/Con)*(R/K/F)
    Result.AR_rt1(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==1,11));
    Result.AK_rt1(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==2,11));
    Result.AF_rt1(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==3 | Abstract_old(:,4)==4,11));
    Result.CR_rt1(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==1,11));
    Result.CK_rt1(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==2,11));
    Result.CF_rt1(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==3 | Concrete_old(:,4)==4,11));

    % Acc: (Abs/Con)*(R/K/F)
    Result.AR_acc(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==1,12));
    Result.AK_acc(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==2,12));
    Result.AF_acc(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==3 | Abstract_old(:,4)==4,12));
    Result.CR_acc(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==1,12));
    Result.CK_acc(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==2,12));
    Result.CF_acc(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==3 | Concrete_old(:,4)==4,12));
    
    % count trial number with different responses
    Result.OldR(s,1) = size(old(old(:,4)==1),1);
    Result.OldK(s,1) = size(old(old(:,4)==2),1);    
    Result.OldU(s,1) = size(old(old(:,4)==3),1);
    Result.OldN(s,1) = size(old(old(:,4)==4),1);
    Result.nOld(s,1) = size(old,1);
    
    Result.NewR(s,1) = size(new(new(:,4)==1),1);
    Result.NewK(s,1) = size(new(new(:,4)==2),1);
    Result.NewU(s,1) = size(new(new(:,4)==3),1);
    Result.NewN(s,1) = size(new(new(:,4)==4),1);
    Result.nNew(s,1) = size(new,1);
    
    Result.AbsR_old(s,1) = size(Abstract_old(Abstract_old(:,4)==1),1);
    Result.AbsK_old(s,1) = size(Abstract_old(Abstract_old(:,4)==2),1);
    Result.AbsU_old(s,1) = size(Abstract_old(Abstract_old(:,4)==3),1);
    Result.AbsN_old(s,1) = size(Abstract_old(Abstract_old(:,4)==4),1);
    Result.nAbs_old(s,1) = size(Abstract_old,1);
    
    Result.ConR_old(s,1) = size(Concrete_old(Concrete_old(:,4)==1),1);
    Result.ConK_old(s,1) = size(Concrete_old(Concrete_old(:,4)==2),1);
    Result.ConU_old(s,1) = size(Concrete_old(Concrete_old(:,4)==3),1);
    Result.ConN_old(s,1) = size(Concrete_old(Concrete_old(:,4)==4),1);
    Result.nCon_old(s,1) = size(Concrete_old,1);
    
    % calculate mean_RT
    Result.RRT_old(s,1) = nanmean(old(old(:,4)==1,6));
    Result.KRT_old(s,1) = nanmean(old(old(:,4)==2,6));
    Result.URT_old(s,1) = nanmean(old(old(:,4)==3,6));
    Result.NRT_old(s,1) = nanmean(old(old(:,4)==4,6));
    
    Result.RRT_new(s,1) = nanmean(new(new(:,4)==1,6));
    Result.KRT_new(s,1) = nanmean(new(new(:,4)==2,6));
    Result.URT_new(s,1) = nanmean(new(new(:,4)==3,6));
    Result.NRT_new(s,1) = nanmean(new(new(:,4)==4,6));
    
    Result.RRT_old_abs(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==1,6));
    Result.KRT_old_abs(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==2,6));
    Result.URT_old_abs(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==3,6));
    Result.NRT_old_abs(s,1) = nanmean(Abstract_old(Abstract_old(:,4)==4,6));
    
    Result.RRT_old_con(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==1,6));
    Result.KRT_old_con(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==2,6));
    Result.URT_old_con(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==3,6));
    Result.NRT_old_con(s,1) = nanmean(Concrete_old(Concrete_old(:,4)==4,6));
    
    % calculate hit/fa/d-primes
    % all images--rem
    Result.R_Hit(s,1) = size(old(old(:,4)==1),1)/size(old,1);
    Result.R_FA(s,1) = size(new(new(:,4)==1),1)/size(new,1);
    if Result.R_Hit(s,1) ~= 0 && Result.R_FA(s,1) ~= 0
        Result.R_dprime(s,1) = norminv(Result.R_Hit(s,1))-norminv(Result.R_FA(s,1));
    elseif Result.R_Hit(s,1) ~= 0 && Result.R_FA(s,1) == 0
        Result.R_dprime(s,1) = norminv(Result.R_Hit(s,1))-norminv(1/(2*size(new,1)));
    else
        Result.R_dprime(s,1) = 0;
    end
    % all images--know
    Result.K_Hit(s,1) = size(old(old(:,4)==2),1)/size(old,1);
    Result.K_FA(s,1) = size(new(new(:,4)==2),1)/size(new,1);
    if Result.K_Hit(s,1) ~= 0 && Result.K_FA(s,1) ~= 0
        Result.K_dprime(s,1) = norminv(Result.K_Hit(s,1))-norminv(Result.K_FA(s,1));
    elseif Result.K_Hit(s,1) ~= 0 && Result.K_FA(s,1) == 0
        Result.K_dprime(s,1) = norminv(Result.K_Hit(s,1))-norminv(1/(2*size(new,1)));
    else
        Result.K_dprime(s,1) = 0;
    end
    % all images--rem+know
    Result.RK_Hit(s,1) = size(old(old(:,4)==1 | old(:,4)==2),1)/size(old,1);
    Result.RK_FA(s,1) = size(new(new(:,4)==1 | new(:,4)==2),1)/size(new,1);
    if Result.RK_Hit(s,1) ~= 0 && Result.RK_FA(s,1) ~= 0
        Result.RK_dprime(s,1) = norminv(Result.RK_Hit(s,1))-norminv(Result.RK_FA(s,1));
    elseif Result.RK_Hit(s,1) ~= 0 && Result.RK_FA(s,1) == 0
        Result.RK_dprime(s,1) = norminv(Result.RK_Hit(s,1))-norminv(1/(2*size(new,1)));
    else
        Result.RK_dprime(s,1) = 0;
    end
    % abstract images
    Result.Abs_Hit(s,1)= size(Abstract_old(Abstract_old(:,4)==1|Abstract_old(:,4)==2),1)/size(Abstract_old,1);
    Result.Abs_FA(s,1)= size(Abstract_new(Abstract_new(:,4)==1|Abstract_new(:,4)==2),1)/size(Abstract_new,1);
    if Result.Abs_Hit(s,1) ~= 0 && Result.Abs_FA(s,1) ~= 0
        Result.Abs_dprime(s,1) = norminv(Result.Abs_Hit(s,1))-norminv(Result.Abs_FA(s,1));
    elseif Result.Abs_Hit(s,1) ~= 0 && Result.Abs_FA(s,1) == 0
        Result.Abs_dprime(s,1) = norminv(Result.Abs_Hit(s,1))-norminv(1/(2*size(new,1)));
    else
        Result.Abs_dprime(s,1) = 0;
    end
    % concrete images
    Result.Con_Hit(s,1)= size(Concrete_old(Concrete_old(:,4)==1|Concrete_old(:,4)==2),1)/size(Concrete_old,1);
    Result.Con_FA(s,1)= size(Concrete_new(Concrete_new(:,4)==1|Concrete_new(:,4)==2),1)/size(Concrete_new,1);
    if Result.Con_Hit(s,1) ~= 0 && Result.Con_FA(s,1) ~= 0
        Result.Con_dprime(s,1) = norminv(Result.Con_Hit(s,1))-norminv(Result.Con_FA(s,1));
    elseif Result.Con_Hit(s,1) ~= 0 && Result.Con_FA(s,1) == 0
        Result.Con_dprime(s,1) = norminv(Result.Con_Hit(s,1))-norminv(1/(2*size(new,1)));
    else
        Result.Con_dprime(s,1) = 0;
    end
  
    clear LSM TSM SM old new Abstract_old Abstract_new Concrete_old Concrete_new
end

save([output_data_path 'Beh_results_AllSub.mat'],'Result')
clear;clc

%% plot variances---hist
clear;clc
cd 'E:\sisi2020\Experiments\Semantic_memory_EEG\ERP_analysis_sisi\Beha_new\beh_results\';
load('Beh_results_AllSub.mat')

% load selected sub list
para_path = 'E:\sisi2020\Experiments\Semantic_memory_EEG\ERP_analysis_sisi\EN_ERP_analysis\ERS_data\within_condition_ERS\ERS_results\';
load([para_path 'selected_sublist.mat']);

SS15_result = Result(index15,:);
SS20_result = Result(index20,:);
SS30_result = Result(index30,:);
% 52,55,58
% plot rem d-prime
figure('Name','Rem d-prime','NumberTitle','off', 'Color','white');
subplot(1,4,1)
histogram(table2array(Result(:,52)));
title('Allsub d-prime')
xlabel('dprime'), ylabel('subject number')

subplot(1,4,2)
histogram(table2array(SS15_result(:,52)));
title('SS15 d-prime')
xlabel('dprime'), ylabel('subject number')

subplot(1,4,3)
histogram(table2array(SS20_result(:,52)));
title('SS20 d-prime')
xlabel('dprime'), ylabel('subject number')

subplot(1,4,4)
histogram(table2array(SS30_result(:,52)));
title('SS30 d-prime')
xlabel('dprime'), ylabel('subject number')
xlabel('dprime'), ylabel('subject number')

% plot know d-prime
figure('Name','Know d-prime','NumberTitle','off', 'Color','white');
subplot(1,4,1)
histogram(table2array(Result(:,55)));
title('Allsub d-prime')
xlabel('dprime'), ylabel('subject number')

subplot(1,4,2)
histogram(table2array(SS15_result(:,55)));
title('SS15 d-prime')
xlabel('dprime'), ylabel('subject number')

subplot(1,4,3)
histogram(table2array(SS20_result(:,55)));
title('SS20 d-prime')
xlabel('dprime'), ylabel('subject number')

subplot(1,4,4)
histogram(table2array(SS30_result(:,55)));
title('SS30 d-prime')
xlabel('dprime'), ylabel('subject number')
xlabel('dprime'), ylabel('subject number')

% plot rem+know d-prime
figure('Name','Rem+Know d-prime','NumberTitle','off', 'Color','white');
subplot(1,4,1)
histogram(table2array(Result(:,58)));
title('Allsub d-prime')
xlabel('dprime'), ylabel('subject number')

subplot(1,4,2)
histogram(table2array(SS15_result(:,58)));
title('SS15 d-prime')
xlabel('dprime'), ylabel('subject number')

subplot(1,4,3)
histogram(table2array(SS20_result(:,58)));
title('SS20 d-prime')
xlabel('dprime'), ylabel('subject number')

subplot(1,4,4)
histogram(table2array(SS30_result(:,58)));
title('SS30 d-prime')
xlabel('dprime'), ylabel('subject number')

output_data_path = 'E:\sisi2020\Experiments\Semantic_memory_EEG\ERP_analysis_sisi\EN_ERP_analysis\ERS_data\within_condition_ERS\ERS_results\select_subs\';
save(fullfile([output_data_path 'Beh_results_plot.mat']))
clear;clc

% % plot RT
% figure('Name','ALLsub RT','NumberTitle','off', 'Color','white');
% cond_list = {'R-RT-old', 'K-RT-old', 'U-RT-old', 'N-RT-old', 'R-RT-new', 'K-RT-new', 'U-RT-new', 'N-RT-new'};
% for i = 1:length(cond_list)
% subplot(2,4,i)
% hist(table2array(Result(:,33+i))');
% title(cond_list{i})
% end


%% plot bar results
% plot mean & SEM of different conditions
clear;clc
cd 'E:\sisi2020\Experiments\Semantic_memory_EEG\ERP_analysis_sisi\EN_ERP_analysis\ERS_data\within_condition_ERS\ERS_results\select_subs\';
load('Beh_results_plot.mat');

% bar(m*n): m--group variance; n--within-group variance
% subgroup(ss15,ss20,ss30,allsub)*d-prime cat(rem,know,rem+know)
SS15_dprime(:,1) = table2array(SS15_result(:,52));
SS15_dprime(:,2) = table2array(SS15_result(:,55));
SS15_dprime(:,3) = table2array(SS15_result(:,58));
SS20_dprime(:,1) = table2array(SS20_result(:,52));
SS20_dprime(:,2) = table2array(SS20_result(:,55));
SS20_dprime(:,3) = table2array(SS20_result(:,58));
SS30_dprime(:,1) = table2array(SS30_result(:,52));
SS30_dprime(:,2) = table2array(SS30_result(:,55));
SS30_dprime(:,3) = table2array(SS30_result(:,58));
AS_dprime(:,1) = table2array(Result(:,52));
AS_dprime(:,2) = table2array(Result(:,55));
AS_dprime(:,3) = table2array(Result(:,58));

% calculate means & SEMs
mean_ss = zeros(4,3);
mean_ss(1,1:3) = mean(SS15_dprime,1);   % use nanmean if NaNs exist
mean_ss(2,1:3) = mean(SS20_dprime,1);
mean_ss(3,1:3) = mean(SS30_dprime,1);
mean_ss(4,1:3) = mean(AS_dprime,1);
% SEM formula: SEM = std(data,2)./sqrt(size(data,2));
SEM_ss = zeros(4,3);
SEM_ss(1,1:3) = std(SS15_dprime,1)./sqrt(size(SS15_dprime,1));
SEM_ss(2,1:3) = std(SS20_dprime,1)./sqrt(size(SS20_dprime,1));
SEM_ss(3,1:3) = std(SS30_dprime,1)./sqrt(size(SS30_dprime,1));
SEM_ss(4,1:3) = std(AS_dprime,1)./sqrt(size(AS_dprime,1));


% plot means with error bars
figure('Name','Beh dprime','NumberTitle','off', 'Color','white'),
hold on
bp = bar(mean_ss);
pause(0.1); %pause allows the figure to be created
% For each set of bars, find the centers of the bars, and write error bars
for ib = 1:numel(bp)
    %XData property is the tick labels/group centers; 
    %XOffset is the offset of each distinct group
    xData = bp(ib).XData+bp(ib).XOffset;
    errorbar(xData,mean_ss(:,ib),SEM_ss(:,ib),'k.')
end
set(gca, 'xtick', 1:size(mean_ss,1), 'xticklabel', {'SS15' 'SS20' 'SS30' 'S231'})
xlabel('Subject group'), ylabel('Dprime')
legend('rem-dprime', 'know-dprime', 'rem+know-dprime')

save(fullfile('Beh_results_plot_bar.mat'));
clear;clc


%%
close all;
clear;clc

%% plot correlation of dprime & ERS left trials
clear;clc
cd 'E:\sisi2020\Experiments\Semantic_memory_EEG\ERP_analysis_sisi\EN_ERP_analysis\ERS_data\within_condition_ERS\ERS_results\select_subs\';
load('Beh_results_plot_bar.mat');
load('ERS_left_trial_stats_AllSub.mat');

x1 = AS_dprime; % Rem, Know, Rem+Know respectively
x2 = SS15_dprime;
x3 = SS20_dprime;
x4 = SS30_dprime;
y1 = left_trial_ERS_WI; % rem, know, for 
y2 = left_trial_sub15;
y3 = left_trial_sub20;
y4 = left_trial_sub30;

%%%%%%%%%%%%% Rem-dprime %%%%%%%%%%%%%%%%%%
figure('Name','Corr Rem-dprime ERS-left-trial','NumberTitle','off', 'Color','white'),
% rem
con_list = {'rem','know','for'};
for i = 1:length(con_list)
subplot(4,3,i)
hold on
scatter(x1(:,1),y1(:,i),[],'r','filled');  % plot dots
xlabel('Rem-dprime'), ylabel('Allsub-left trial number')
title(con_list{i})
end

for i = 1:length(con_list)
subplot(4,3,i+3)
hold on
scatter(x2(:,1),y2(:,i),[],'m','filled');  % plot dots
xlabel('Rem-dprime'), ylabel('SS15-left trial number')
title(con_list{i})
end

for i = 1:length(con_list)
subplot(4,3,i+6)
hold on
scatter(x3(:,1),y3(:,i),[],'b','filled');  % plot dots
xlabel('Rem-dprime'), ylabel('SS20-left trial number')
title(con_list{i})
end

for i = 1:length(con_list)
subplot(4,3,i+9)
hold on
scatter(x4(:,1),y4(:,i),[],'k','filled');  % plot dots
xlabel('Rem-dprime'), ylabel('SS30-left trial number')
title(con_list{i})
end

%%%%%%%%%%%%% Know-dprime %%%%%%%%%%%%%%%%%%
figure('Name','Corr Know-dprime ERS-left-trial','NumberTitle','off', 'Color','white'),
% rem
con_list = {'rem','know','for'};
for i = 1:length(con_list)
subplot(4,3,i)
hold on
scatter(x1(:,2),y1(:,i),[],'r','filled');  % plot dots
xlabel('Know-dprime'), ylabel('Allsub-left trial number')
title(con_list{i})
end

for i = 1:length(con_list)
subplot(4,3,i+3)
hold on
scatter(x2(:,2),y2(:,i),[],'m','filled');  % plot dots
xlabel('Know-dprime'), ylabel('SS15-left trial number')
title(con_list{i})
end

for i = 1:length(con_list)
subplot(4,3,i+6)
hold on
scatter(x3(:,2),y3(:,i),[],'b','filled');  % plot dots
xlabel('Know-dprime'), ylabel('SS20-left trial number')
title(con_list{i})
end

for i = 1:length(con_list)
subplot(4,3,i+9)
hold on
scatter(x4(:,2),y4(:,i),[],'k','filled');  % plot dots
xlabel('Know-dprime'), ylabel('SS30-left trial number')
title(con_list{i})
end

%%%%%%%%%%%%% Rem+Know-dprime %%%%%%%%%%%%%%%%%%
figure('Name','Corr Rem+Know-dprime ERS-left-trial','NumberTitle','off', 'Color','white'),
% rem
con_list = {'rem','know','for'};
for i = 1:length(con_list)
subplot(4,3,i)
hold on
scatter(x1(:,3),y1(:,i),[],'r','filled');  % plot dots
xlabel('RemKnow-dprime'), ylabel('Allsub-left trial number')
title(con_list{i})
end

for i = 1:length(con_list)
subplot(4,3,i+3)
hold on
scatter(x2(:,3),y2(:,i),[],'m','filled');  % plot dots
xlabel('RemKnow-dprime'), ylabel('SS15-left trial number')
title(con_list{i})
end

for i = 1:length(con_list)
subplot(4,3,i+6)
hold on
scatter(x3(:,3),y3(:,i),[],'b','filled');  % plot dots
xlabel('RemKnow-dprime'), ylabel('SS20-left trial number')
title(con_list{i})
end

for i = 1:length(con_list)
subplot(4,3,i+9)
hold on
scatter(x4(:,3),y4(:,i),[],'k','filled');  % plot dots
xlabel('RemKnow-dprime'), ylabel('SS30-left trial number')
title(con_list{i})
end

hold off
clear;clc

%% select subs with high Know-dprime
clear;clc
cd 'E:\sisi2020\Experiments\Semantic_memory_EEG\ERP_analysis_sisi\EN_ERP_analysis\ERS_data\within_condition_ERS\ERS_results\select_subs\';
load('Beh_results_plot_bar.mat');
load('ERS_left_trial_stats_AllSub.mat');

index_Kdprime.SS15 = index15(find(SS15_dprime(:,2)>0.5));
index_Kdprime.SS20 = index20(find(SS20_dprime(:,2)>0.5));
index_Kdprime.SS30 = index30(find(SS30_dprime(:,2)>0.5));

save('SS_index_Kdprime.mat', 'index_Kdprime')
clear;clc

