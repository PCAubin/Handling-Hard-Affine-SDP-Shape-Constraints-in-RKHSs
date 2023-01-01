function [X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq,...
    D_arr_uniq,D_cellArr_uniq_Ineq,bool_arr_C_Ineq_within_C_sol] =...
    ConvertingConsProblem_toUniquePoints_wDeriv(X_arr,Y_arr,C_arr,D_arr,X_grid,Y_grid)
% ConvertingConsProblem_toUniquePoints transforms the inputs
% defining the constrained problem as given within the function
% "..._generation" into an output for which the evaluation points are
% unique (to have a better conditionned matrix [K(x_i,x_j)]_{i,j}). The
% inputs X_arr,Y_arr,C_arr,D_arr are cell arrays 3x1, X_grid,Y_grid are matrices.
% For the input the convention is:
% 1: objective points, 2: inequality points, 3: equality points
% For the output the convention is:
% 1: solution points, 2: objective points, 3: inequality points,
% 4: equality points, 5: grid points


Xtot=cell2mat(X_arr); Ctot=cell2mat(C_arr); Ytot=cell2mat(Y_arr); Dtot=cell2mat(D_arr);
n_arr=[0;size(X_arr{1},1);...
    size(X_arr{2},1)+size(X_arr{1},1);...
    size(X_arr{3},1)+size(X_arr{2},1)+size(X_arr{1},1);size(Xtot,1)];

if ~(size(Ytot,1)==size(Xtot,1))
   error(['Dimension mismatch occurred: Total number of points is ', ...
            int2str(size(Xtot,1)),' while total number of values is ', ...
            int2str(size(Ytot,1)),'.'])
end
if ~(size(Ctot,1)==size(Xtot,1))
   error(['Dimension mismatch occurred: Total number of points is ', ...
            int2str(size(Xtot,1)),' while total number of constraints is ', ...
            int2str(size(Ctot,1)),'.'])
end

[Xtot_uniq,~,idx_XuniqXtot]=unique(Xtot,'rows','stable'); Ntot_uniq = length(Xtot_uniq);

sub_YInUniq=cell(Ntot_uniq,4);
sub_CInUniq=cell(Ntot_uniq,4);
sub_DInUniq=cell(Ntot_uniq,4);

X_arr_uniq=cell(5,1);
Y_arr_uniq=cell(5,1);
C_arr_uniq=cell(4,1);
D_arr_uniq=cell(4,1);

for i = 1:Ntot_uniq
      sub_CInUniq{i,1} = sparse(Ctot(idx_XuniqXtot==i,:)');
      sub_DInUniq{i,1} = sparse(Dtot(idx_XuniqXtot==i,:)');
      sub_YInUniq{i,1} = sparse(Ytot(idx_XuniqXtot==i,:)');
%        [~,pivot_col] = rref(sub_CInUniq{i,1});%IN CASE THE CONSTRAINTS AT POINT Xi ARE NOT INDEPENDENT
%       sub_CInUniq{i,1}=sub_CInUniq{i,1}(:,pivot_col);
%       sub_DInUniq{i,1}=sub_DInUniq{i,1}(:,pivot_col);
end

C_arr_uniq{1}= blkdiag(sub_CInUniq{:,1});
% D_arr_uniq{1}= blkdiag(sub_DInUniq{:,1});
D_arr_uniq{1}= cell2mat(sub_DInUniq(:,1)');
Y_arr_uniq{1}= blkdiag(sub_YInUniq{:,1});
X_arr_uniq{1}= Xtot_uniq;    

X_arr_uniq{end}= X_grid; 
Y_arr_uniq{end}= Y_grid; 

for j=1:3
    for i = 1:Ntot_uniq
        temp_idx=(idx_XuniqXtot==i)&((1:n_arr(end))'>n_arr(j))...
            &((1:n_arr(end))'<=n_arr(j+1));
      sub_CInUniq{i,j+1} = sparse(Ctot(temp_idx,:)');
      sub_DInUniq{i,j+1} = sparse(Dtot(temp_idx,:)');
      sub_YInUniq{i,j+1} = sparse(Ytot(temp_idx,:)');
    end
    X_arr_uniq{j+1}= Xtot_uniq(ismember(Xtot_uniq,X_arr{j},'rows'),:);
    Y_arr_uniq{j+1}= blkdiag(sub_YInUniq{:,j+1});
    C_arr_uniq{j+1}= blkdiag(sub_CInUniq{:,j+1});
    D_arr_uniq{j+1}= cell2mat(sub_DInUniq(:,j+1)');
end
C_cellArr_uniq_Ineq= sub_CInUniq(:,3);
C_cellArr_uniq_Ineq=C_cellArr_uniq_Ineq(~cellfun('isempty',C_cellArr_uniq_Ineq));
D_cellArr_uniq_Ineq= sub_DInUniq(:,3);
D_cellArr_uniq_Ineq=D_cellArr_uniq_Ineq(~cellfun('isempty',D_cellArr_uniq_Ineq));

%ATTEMPT AT RETRIEVING THE INDICES BUT JORDAN DECOMPOSITION WOULD MESS IT
%UP
% idx_arr_C_Ineq_within_C_sol=false(size(C_arr_uniq{1},2),1);
sub_CInUniq_idx=cell(Ntot_uniq,1);
sub_CInUniq_Ineq_idx=cell(Ntot_uniq,1);
for i = 1:Ntot_uniq
    j=2;
    sub_CInUniq_idx{i,1} = find(idx_XuniqXtot==i);
    sub_temp_idx=find((idx_XuniqXtot==i)&((1:n_arr(end))'>n_arr(j))...
            &((1:n_arr(end))'<=n_arr(j+1)));
    sub_CInUniq_Ineq_idx{i,1}=sub_temp_idx;
        %[~,idx]=ismember(B,A)
%      idx_arr_C_Ineq_within_C_sol(count:count+size(temp_idx)) = sparse(Ctot(idx_XuniqXtot==i,:)');
%      count=
end
[~,idx_arr_C_Ineq_within_C_sol]=ismember(cell2mat(sub_CInUniq_Ineq_idx),cell2mat(sub_CInUniq_idx));
bool_arr_C_Ineq_within_C_sol=false(size(C_arr_uniq{1},2),1);
bool_arr_C_Ineq_within_C_sol(idx_arr_C_Ineq_within_C_sol)=true;
end

