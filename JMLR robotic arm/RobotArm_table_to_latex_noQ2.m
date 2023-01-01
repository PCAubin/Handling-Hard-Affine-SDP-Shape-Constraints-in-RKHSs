% We use latexTable.m to generate automatically a latex table of the results stored in
% table_results_RobotArm_4D_fullGaussian_wErrBdary_500ptEval.txt. 

mean_matErr=[]; std_matErr=[]; rowIdx=[];
for dim_in=[4 6]
err_arr = dlmread(['table_results_RobotArm_' num2str(dim_in) 'D_fullGaussian_wErrBdary_500ptEval.txt'],'\t', 0,0);
error_list=[{'L_{\text{cons}}^2'};{'L_{\mathrm{cons}}^1'};{'T_s'}];
error_idx_list=[3 5 6];
error_idx_wTime_list=[3 6 17];

N_array=unique(err_arr(:,7));%5:5:15;
DFac_array=unique(err_arr(:,9));

err_arr=err_arr(err_arr(:,9)==50,:);
%%ONE BIG MATRIX (grouped by discretization)
error_list_temp=error_list;
mean_mat=nan(7,length(N_array)*length(error_list_temp));
std_mat=nan(7,length(N_array)*length(error_list_temp));
for idxN=1:length(N_array)
    N=N_array(idxN);
    mean_mat(1,1+(idxN-1)*length(error_list_temp):idxN*length(error_list_temp))=N^dim_in;
    err_arr_temp=err_arr(err_arr(:,7)==N,:);
    mean_mat(2,1+(idxN-1)*length(error_list_temp):idxN*length(error_list_temp))=ceil(nanmean(err_arr_temp(:,8)));
for idxMeth=0:4
    sub_err_arr=err_arr_temp(err_arr_temp(:,1)==idxMeth,error_idx_wTime_list);
    mean_mat(3+idxMeth,1+(idxN-1)*length(error_list_temp):idxN*length(error_list_temp))=nanmean(sub_err_arr);
    std_mat(3+idxMeth,1+(idxN-1)*length(error_list_temp):idxN*length(error_list_temp))=nanstd(sub_err_arr);
%     std_mat(3+idxMeth,idxN*length(error_list_temp))=0;
end
end
mean_mat=[[dim_in,nan(1,size(mean_mat,2)-1)];mean_mat];
std_mat=[nan(1,size(mean_mat,2));std_mat];
% mean_mat=table(error_listmean_mat);
%PERMUTE COLUMNS
reIdx=repmat([[0:length(N_array)-1]*length(error_list_temp)]',1,length(error_list_temp))+repmat(1:length(error_list_temp),length(N_array),1);
mean_matErr=[mean_matErr;mean_mat(:,reshape(reIdx,1,[]))'];
std_matErr=[std_matErr;std_mat(:,reshape(reIdx,1,[]))'];
rowIdx=[rowIdx;length(N_array)];
end
%% GENERATE LATEX TABLE
rowNames=cell(size(mean_matErr,1),1);
rowNames(:) = {''};
shortRowNames=[{'$L_{\text{cons}}^2$'};{'$L_{\mathrm{cons}}^1$'};{'$T_s$'}];
startIdx=1;
for count=1:length(rowIdx)
    for i=1:length(error_list_temp)
        rowNames{startIdx+(i-1)*rowIdx(count)}=shortRowNames{i};
    end
    startIdx=startIdx+i*rowIdx(count);
end
columnNames={'\begin{tabular}{c} $d$ \\ dim\end{tabular}';...
    '\begin{tabular}{c}$M$\\ \# pts\end{tabular}';...
    '\begin{tabular}{c}$\approx d\times M$\\ \# cons\end{tabular}';...
    'Unconstrained';'\begin{tabular}{c}Discretized\\ constraints\end{tabular}';...
    '\begin{tabular}{c}SOC-ball\\ constraints\end{tabular}';...
    '\begin{tabular}{c}SOC-hyp\\ constraints\end{tabular}';...
    '\begin{tabular}{c}kSoS\\ constraints\end{tabular}'};
% Clear the selected options from previous example
clear input;
fprintf('\n\nUsing string data that includes LaTex code in a MATLAB table\n\n');
complete_pm_table=cell(size(mean_matErr));
for i=1:size(mean_matErr,1)
for j=1:size(mean_matErr,2)
    if j<4
    complete_pm_table{i,j}=num2str(mean_matErr(i,j),'%d');
    else
        if mean_matErr(i,j)==0
            complete_pm_table{i,j}='$<$0.01';
        elseif isnan(mean_matErr(i,j))
            complete_pm_table{i,j}='--';
        elseif isnan(std_matErr(i,j))||(std_matErr(i,j)==0)
            complete_pm_table{i,j}=sprintf('%.3f',mean_matErr(i,j));
        else
            complete_pm_table{i,j}=[sprintf('%.3f',mean_matErr(i,j)) ' $\pm$ ' sprintf('%1.0e',std_matErr(i,j))];
        end
    end
end     
end
input.data=table([columnNames';complete_pm_table]);
% Set column labels (use empty string for no label):
% input.tableColLabels = columnNames';
% Set row labels (use empty string for no label):
input.tableRowLabels = [{''};rowNames]';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;
% input.data=mean_matErr;
input.dataNanString = '--';
input.tableColumnAlignment = 'l';
input.tableBorders = 1;
input.booktabs = 1;
% Now call the function to generate LaTex code:
latex = latexTable(input);
fid=fopen('MyLatex_pm_noQ2.tex','w');
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
fprintf('\n... your LaTex code has been saved as ''MyLatex_pm_noQ2.tex'' in your working directory\n');

%{' & dim & \# pts & \# cons & Unconstrained & Sampled constraints & SOC-ball constraints & SOC-hyp constraints & kSoS constraints \\'}
