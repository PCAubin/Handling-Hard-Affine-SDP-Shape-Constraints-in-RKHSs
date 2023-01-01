% The code below outputs the boxplots of the article based on the
% table of results of the economy experiment.


MSE_arr = dlmread('KRR_Labour_MSEVals_20splits_4methods_final.txt','\t', 0,0);
figure
hold on
color_list=['r';'m';'b'];
list_tags=["NoCons"; "SOC Conv."; "SOC Monot." ; "SOC Conv.+Monot."];
boxplot(reshape(MSE_arr(:,1),20,[]),'Color','r', 'Widths',0.2)
boxplot(reshape(MSE_arr(:,3),20,[]),'Color','b', 'Widths',0.2)%,'Labels',list_tags
set(gca,'XTickLabel',list_tags,'FontSize',34)

box_vars = findobj(gca,'Tag','Box');
legend(box_vars([1,5]),[{'Train (27 points - 5% of total)'},{'Test (272 points - 50% of total)'}],...
    'Interpreter','latex','FontSize',34);
ylabel('Mean-Squared Error (MSE)','FontSize',34)   
ylim([0 0.65])

% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,['KRR_LabourData_Boxplots_4methods_large'],'-dpdf')

%% Small plot
figure
hold on
list_tags=["NoCons"; "SOC C"; "SOC M" ; "SOC C+M"];
boxplot(reshape(MSE_arr(:,1),20,[]),'Color','r')
boxplot(reshape(MSE_arr(:,3),20,[]),'Color','b')%,'Labels',list_tags
set(gca,'XTickLabel',list_tags,'FontSize',12)

box_vars = findobj(gca,'Tag','Box');
legend(box_vars([1,5]),[{'Train (27 points - 5% of total)'},{'Test (272 points - 50% of total)'}],...
    'Interpreter','latex','FontSize',14);
ylabel('Mean-Squared Error (MSE)','FontSize',12) 
ylim([0 0.65])
% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,['KRR_LabourData_Boxplots_4methods_small'],'-dpdf')