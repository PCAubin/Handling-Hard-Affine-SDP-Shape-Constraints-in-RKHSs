% The code below allows to plot some figures to study the content of the
% table of results for the robot arm experiment.

err_arr = dlmread('table_results_RobotArm_6D_fullGaussian_wErrBdary_500ptEval.txt','\t', 0,0);
color_list=['g';'r';'b';'c';'m';'k'];
% lambdaSoS_list=1E-3; %[1E-6;1E-5;1E-4;5E-4;1E-3;5E-3;1E-2];
error_list=[{'Q2 error'};{'L2 error global'};{'L_{infty} error constraints'};{'L_1 error constraints'}];
error_idx_list=[2 3 5 6];
legend_entries=[{'No constraints'},{'Sampled constraints'},{'SOCball constraints'},...
    {'SOChyp constraints'},{'kSoS constraints'},{'Zero solution'}];

% Ndisc=2;
% err_arr=err_arr(err_arr(:,7)==Ndisc,:);
N_array=unique(err_arr(:,7));%5:5:15;
%
ax_bounds=[min(err_arr(:,2)) max(err_arr(:,2));...
    min(err_arr(:,3)) max(err_arr(:,3));...
    min(err_arr(:,4)) max(err_arr(:,4));...
%     min(err_arr(:,6)) max(err_arr(:,6))];
    log10(min(err_arr(:,6))) log10(max(err_arr(:,6)))];

ax_bounds=[-0.1 1.1;...
    0 1.5;...
    0 3;...
    min(err_arr(:,6)) max(err_arr(:,6))];
%     log10(min(err_arr(:,6))) log10(max(err_arr(:,6)))];


close all
% idxLSoS=4;
% lambdaSoS=lambdaSoS_list(idxLSoS);


for idxN=1:length(N_array)
    figure(idxN)
for idxE=1:length(error_list)
%1:length(N_array)%5:10:15
    N=N_array(idxN);
%     lambdaSoS=lambdaSoS_list(idxLSoS);
err_arr_temp=err_arr(err_arr(:,7)==N,:);
M_array=sort(unique(err_arr_temp(:,9)))';

subplot(1,length(error_list),idxE)
count_error=error_idx_list(idxE);
% count_color=1;
title([error_list(idxE)])%'N=' num2str(N), 
hold on
idx_array=sort(unique(err_arr_temp(:,1)));
for j=1:length(idx_array)
sub_idx=(err_arr_temp(:,1)==idx_array(j));
sub_mat=err_arr_temp(sub_idx,:);
mat_to_plot=nan(size(sub_mat,1),length(M_array));
mat_to_plot(any(isnan(mat_to_plot),2),:) = [];
for i=1:length(M_array)
    data=sub_mat(sub_mat(:,9)==M_array(i),count_error);
    mat_to_plot(1:size(data,1),i)=data+1E-20;    
end    

% if (error_idx_list(idxE)==6)%||(error_idx_list(idxE)==3)
% %     sum(find(mat_to_plot==0))
% %     mat_to_plot(mat_to_plot==0)=NaN;
%     mat_to_plot=log10(mat_to_plot+1E-20);%
% end
line(1:(length(M_array)),nanmedian(mat_to_plot),'Color',color_list(idx_array(j)+1),'LineWidth',2 )
boxplot(mat_to_plot,'Labels',...
    {num2str(M_array')},'Color',color_list(idx_array(j)+1))
% count_color=count_color+1;
end
xlabel('M')
ylabel(error_list(count_error==error_idx_list))
legend(legend_entries(idx_array+1),...
    'Interpreter','latex');%,'FontSize',34
ylim(ax_bounds(idxE,:))
end
end
