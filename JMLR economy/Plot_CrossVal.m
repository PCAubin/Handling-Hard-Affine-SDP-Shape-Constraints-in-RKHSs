Ploss_arr1 = dlmread('CrossVal_table_results_Labour.txt','\t', 0,0);
Ploss_arr2 = dlmread('CrossVal_table_results_Labour_wNorm.txt','\t', 0,0);
close all
% figure
% hold on
% idx=(Ploss_arr1(:,end-1)==0);
% log_arr=log10(Ploss_arr1(idx,7));
% plot(log_arr,Ploss_arr1(idx,3),'+g')%train
% plot(log_arr,Ploss_arr1(idx,1),'+b')%test
% plot(log_arr,Ploss_arr1(idx,2),'+r')%testcomp
% 
% idx=(Ploss_arr1(:,end-1)==1);
% log_arr=log10(Ploss_arr1(idx,7));
% plot(log_arr,Ploss_arr1(idx,3),'og')
% plot(log_arr,Ploss_arr1(idx,1),'ob')
% plot(log_arr,Ploss_arr1(idx,2),'or')
% 
% log_arr2=log10(Ploss_arr2(:,8));
% plot(log_arr2,Ploss_arr2(:,3),'*g')
% plot(log_arr2,Ploss_arr2(:,1),'*b')
% plot(log_arr2,Ploss_arr2(:,2),'*r')
% ylim([0,3])
%% No constraints
figure
hold on

idx=(Ploss_arr1(:,end-1)==0);
L_array=sort(unique(Ploss_arr1(idx,7)))';

mat_to_plot=nan(size(Ploss_arr1,1),length(L_array));
mat_to_plot(any(isnan(mat_to_plot),2),:) = [];
color_list=['b';'r';'m'];
for j=[3,2,1]
for i=1:length(L_array)
    tmpData=Ploss_arr1(Ploss_arr1(:,7)==L_array(i),j);
    mat_to_plot(1:size(tmpData,1),i)=tmpData;    
end    
mat_to_plot(mat_to_plot==0)=NaN;
line(1:(length(L_array)),nanmean(mat_to_plot),'Color',color_list(j),'LineWidth',2)
boxplot(mat_to_plot,'Color',color_list(j),'Labels',...
    {num2str(L_array')})
end
ylim([0 1])
title('No constraints')
%% All discretized constraints
figure
hold on

idx=(Ploss_arr1(:,end-1)==1);
L_array=sort(unique(Ploss_arr1(idx,7)))';

mat_to_plot=nan(size(Ploss_arr1,1),length(L_array));
mat_to_plot(any(isnan(mat_to_plot),2),:) = [];
color_list=['b';'r';'m'];
for j=[3,2,1]
for i=1:length(L_array)
    tmpData=Ploss_arr1(Ploss_arr1(:,7)==L_array(i),j);
    mat_to_plot(1:size(tmpData,1),i)=tmpData;    
end    
mat_to_plot(mat_to_plot==0)=NaN;
line(1:(length(L_array)),nanmean(mat_to_plot),'Color',color_list(j),'LineWidth',2)
boxplot(mat_to_plot,'Color',color_list(j),'Labels',...
    {num2str(L_array')})
end
title('Both discretized constraints')
ylim([0 1])
%% Discretized monotonicity
figure
hold on

idx=all(Ploss_arr2(:,end-3:end)==[1  	 0  	 0  	 0 ],2);
L_array=sort(unique(Ploss_arr2(idx,8)))';
mat_to_plot=nan(size(Ploss_arr2,1),length(L_array));
mat_to_plot(any(isnan(mat_to_plot),2),:) = [];
color_list=['b';'r';'m'];
for j=[3,2,1]
for i=1:length(L_array)
    tmpData=Ploss_arr2(Ploss_arr2(:,8)==L_array(i),j);
    mat_to_plot(1:size(tmpData,1),i)=tmpData;    
end    
mat_to_plot(mat_to_plot==0)=NaN;
line(1:(length(L_array)),nanmean(mat_to_plot),'Color',color_list(j),'LineWidth',2)
boxplot(mat_to_plot,'Color',color_list(j),'Labels',...
    {num2str(L_array')})
end
ylim([0 1])
title('Discretized monotonicity')
xlabel('lambda')
ylabel('val on test')
%% SOC monotonicity
figure
hold on

idx=all(Ploss_arr2(:,end-3:end)==[1  	 1  	 0  	 0 ],2);
L_array=sort(unique(Ploss_arr2(idx,8)))';
mat_to_plot=nan(size(Ploss_arr2,1),length(L_array));
mat_to_plot(any(isnan(mat_to_plot),2),:) = [];
color_list=['b';'r';'m'];
for j=[3,2,1]
for i=1:length(L_array)
    tmpData=Ploss_arr2(Ploss_arr2(:,8)==L_array(i),j);
    mat_to_plot(1:size(tmpData,1),i)=tmpData;    
end    
mat_to_plot(mat_to_plot==0)=NaN;
line(1:(length(L_array)),nanmean(mat_to_plot),'Color',color_list(j),'LineWidth',2)
boxplot(mat_to_plot,'Color',color_list(j),'Labels',...
    {num2str(L_array')})
end
ylim([0 1])
title('SOC monotonicity')
xlabel('lambda')
ylabel('val on test')

%% Both SOC constraints
figure
hold on

idx=all(Ploss_arr2(:,end-3:end)==[1  	 1  	 1 	 1 ],2);
L_array=sort(unique(Ploss_arr2(idx,8)))';
mat_to_plot=nan(size(Ploss_arr2,1),length(L_array));
mat_to_plot(any(isnan(mat_to_plot),2),:) = [];
color_list=['b';'r';'m'];
for j=[3,2,1]
for i=1:length(L_array)
    tmpData=Ploss_arr2(Ploss_arr2(:,8)==L_array(i),j);
    mat_to_plot(1:size(tmpData,1),i)=tmpData;    
end    
mat_to_plot(mat_to_plot==0)=NaN;
line(1:(length(L_array)),nanmean(mat_to_plot),'Color',color_list(j),'LineWidth',2)
boxplot(mat_to_plot,'Color',color_list(j),'Labels',...
    {num2str(L_array')})
end
ylim([0 1])
title('Both SOC constraints')
xlabel('lambda')
ylabel('val on test')