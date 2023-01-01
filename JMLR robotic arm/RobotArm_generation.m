function [X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq,...
    D_arr_uniq,D_cellArr_uniq_Ineq,X_consIneq_eval,C_consIneq_eval,...
    D_consIneq_eval,bool_arr_C_Ineq_within_C_sol] =...
    RobotArm_generation(n_grid,n_samples,n_consIneq_eval,n_consIneq_disc,sigma_noise,n_segments,boolGrid_consIneq_disc)
% RobotArm_generation takes the number of points of the various quantities
% appearing in the problem and outputs the generated points and constraint
% vectors.
% THE PARAMETERS OF THE CONTROL PROBLEM SHOULD BE SET HERE MANUALLY. To
% each unique point x_m (wherever it appears in the problem) corresponds
% a list of c_{j,m}^\top D_{j,m} f(x_m) appearing in the problem. For each c_{j,m} there
% is a real value y_{j,m}, e.g. in the objective y_{j,m}-c_{j,m}^\top D_j f(x_m) or
% in the constraints y_{j,m}+c_{j,m}^\top D_j f(x_m) \ge || = 0. The points
% should be written with repetition, and set in a matrix "X_..". To each
% corresponds a vector c_{j,m} Nx1, stored as a row of "C_..", a real y_{j,m}
% stored as a row of "Y_..", and d_{j,m} the index of the derivatives
% involved, stored as a row of "D_.."
%  The convention for the order in "..._arr" is:
% 1: objective points, 2: inequality points, 3: equality points

% For the output "..._arr_unique" the convention is:
% 1: solution points, 2: objective points, 3: inequality points,
% 4: equality points, 5: grid points

% In the output, the evaluation points are
% unique (to have a better conditionned matrix [K(x_i,x_j)]_{i,j}). The
% inputs "X,Y_arr_unique" are cell arrays 5x1, "C and D_arr_unique" are 4x1. 
% Each cell of "Y,C_arr_unique" contains a block diagonal matrix of the
% reals/vectors to compute y_{j,m}+c_{j,m}^\top D_j f(x_m) when calling the function
% "..._solver" (WHERE OTHER PARAMETERS SHOULD BE SET)

% C_cellArr_uniq_Ineq is a cell_array containing the c_{j,m} for each unique x_m of the
% inequality constraints applied to it. 

% noSmallVals_bool=true removes the points for which c_{j,m} is small.
% These lead to numerical problems or the absence of non-trivial solution, by forcing
% the function to be both nonnegative and nonpositive someewhere.

% n_grid=10;n_samples=20;n_consIneq_eval=200;n_consIneq_disc=100;sigma_noise=0;n_segments=2;
dim_in=2*n_segments;
% if mod(dim_in,2)~=0
%     msg = 'Input dimension is not even';
%     error(msg)
% end
dim_out=3;
noSmallVals_bool=true;

n_constraint_types=2; bool_obj=true; bool_grid=true;

X_arr=cell(bool_obj+bool_grid+n_constraint_types,1);
Y_arr=cell(bool_obj+bool_grid+n_constraint_types,1);
C_arr=cell(bool_obj+n_constraint_types,1);
D_arr=cell(bool_obj+n_constraint_types,1);

% GENERATES SAMPLES FROM A TRAJECTORY
X_samples=lhsdesign(n_samples,dim_in,'Criterion','maximin','Iterations',100);
Y_samples=robot_arm(X_samples);
Y_samples=Y_samples+sigma_noise*randn(n_samples,dim_out);
%[y_design,z_design]
% cov(normalize([y_design,z_design]))

Y_obj=Y_samples(:);%[y_design(:,1);y_design(:,2);z_design];
X_obj=repmat(X_samples,dim_out,1);
C_obj=[[ones(n_samples,1),zeros(n_samples,1),zeros(n_samples,1)];...
    [zeros(n_samples,1),ones(n_samples,1),zeros(n_samples,1)];...
    [zeros(n_samples,1),zeros(n_samples,1),ones(n_samples,1)]];
D_obj=sparse(size(C_obj,1),dim_in);

if boolGrid_consIneq_disc
    nbPts_per_dimCons=max(2,floor(nthroot(n_consIneq_disc,dim_in)));
    deltaDesi=1/nbPts_per_dimCons;
    xCoord=deltaDesi*(1/2+(0:nbPts_per_dimCons-1));
    stringCommandLeft='['; stringCommandRight2='[';
    for i=1:dim_in
        stringCommandLeft=[stringCommandLeft,['x_ineq',num2str(i),',']];
        stringCommandRight2=[stringCommandRight2,['x_ineq',num2str(i),'(:),']];
    end
    stringCommandLeft=[stringCommandLeft(1:end-1),']'];
    stringCommandRight2=[stringCommandRight2(1:end-1),']'];
    stringCommandRight=['= ndgrid(',repmat('xCoord,',1,n_segments),repmat('xCoord,',1,n_segments)]; %2*pi*
    stringCommandRight=[stringCommandRight(1:end-1),');'];
    eval([stringCommandLeft,stringCommandRight]);
    X_consIneq_orig=eval(stringCommandRight2);
else
X_consIneq_orig=lhsdesign(n_consIneq_disc,dim_in,'Criterion','maximin','Iterations',100);
end
X_consIneq_orig=X_consIneq_orig+(rand(size(X_consIneq_orig))-ones(size(X_consIneq_orig))/2)*1E-5;
theta_sums = 2*pi*cumsum(X_consIneq_orig(:,n_segments+1:end),2);
X_consIneq=reshape(repmat(X_consIneq_orig',dim_in,1),dim_in,[])';
n_consIneq_tot=size(X_consIneq,1);
C_consIneq=zeros(n_consIneq_tot,dim_out);
D_consIneq=zeros(n_consIneq_tot,dim_in);
for i=1:size(X_consIneq_orig,1)
    C_consIneq((i-1)*dim_in+1:i*dim_in-n_segments,1)=cos(theta_sums(i,1:n_segments))';
    C_consIneq((i-1)*dim_in+1+n_segments:i*dim_in,2)=sin(theta_sums(i,1:n_segments))';
	% ALTERNATIVE WHERE ONLY THE SIGN IS CONSIDERED, IN PRACTICE THIS CAUSED SOME NUMERICAL INSTABILITY WHEN SOLVING
%     C_consIneq((i-1)*dim_in+1:i*dim_in-n_segments,1)=sign(cos(theta_sums(i,1:n_segments))');
%     C_consIneq((i-1)*dim_in+1+n_segments:i*dim_in,2)=sign(sin(theta_sums(i,1:n_segments))');
    D_consIneq((i-1)*dim_in+1:i*dim_in-n_segments,1:n_segments)=eye(n_segments);
    D_consIneq((i-1)*dim_in+1+n_segments:i*dim_in,1:n_segments)=eye(n_segments);
end
% C_consIneq=sparse(C_consIneq);

% CORRESPONDING REAL VALUES AT CONSTRAINT (d_i(x_m))
Y_consIneq=sparse(n_consIneq_tot,1);

%REMOVE SMALL VALUES
if noSmallVals_bool
threshold_Vals=1E-1;
[rows_to_remove, ~] = find((abs(C_consIneq)<threshold_Vals)&(C_consIneq~=0));
rows_to_remove=unique(rows_to_remove);
disp(['Number of rows removed to avoid too small constraint values: '...
    num2str(length(rows_to_remove)) ' over ' num2str(size(C_consIneq,1)) ' rows'])
C_consIneq(rows_to_remove,:)=[];
D_consIneq(rows_to_remove,:)=[];
X_consIneq(rows_to_remove,:)=[];
Y_consIneq(rows_to_remove,:)=[];
end

% X_consIneq=unique(X_consIneq,'rows');

% n_consIneq_eval=10;n_consIneq_disc=5;
%D_arr_uniq,D_consIneq_eval
%Evaluation points for violation of constraints
X_consIneq_eval_orig=lhsdesign(n_consIneq_eval,dim_in,'Criterion','maximin','Iterations',100);
theta_sums = 2*pi*cumsum(X_consIneq_eval_orig(:,n_segments+1:end),2);
X_consIneq_eval=reshape(repmat(X_consIneq_eval_orig',dim_in,1),dim_in,[])';
n_consIneq_eval_tot=size(X_consIneq_eval,1);
C_consIneq_eval=zeros(n_consIneq_eval_tot,dim_out);
D_consIneq_eval=zeros(n_consIneq_eval_tot,dim_in);
for i=1:size(X_consIneq_eval_orig,1)
    C_consIneq_eval((i-1)*dim_in+1:i*dim_in-n_segments,1)=cos(theta_sums(i,1:n_segments))';
    C_consIneq_eval((i-1)*dim_in+1+n_segments:i*dim_in,2)=sin(theta_sums(i,1:n_segments))';
%     C_consIneq_eval((i-1)*dim_in+1:i*dim_in-n_segments,1)=sign(cos(theta_sums(i,1:n_segments))');
%     C_consIneq_eval((i-1)*dim_in+1+n_segments:i*dim_in,2)=sign(sin(theta_sums(i,1:n_segments))');  
    D_consIneq_eval((i-1)*dim_in+1:i*dim_in-n_segments,1:n_segments)=eye(n_segments);
    D_consIneq_eval((i-1)*dim_in+1+n_segments:i*dim_in,1:n_segments)=eye(n_segments);
end
% ZERO VECTOR AT ORIGIN
X_consZ=zeros(dim_out,dim_in);
Y_consZ=zeros(dim_out,1);
C_consZ=eye(dim_out);
D_consZ=sparse(size(C_consZ,1),dim_in);

X_arr{1}=X_obj; Y_arr{1}=Y_obj; C_arr{1}=C_obj; D_arr{1}=D_obj;
X_arr{2}=X_consIneq; Y_arr{2}=Y_consIneq; C_arr{2}=C_consIneq; D_arr{2}=D_consIneq;% Inequality C^\top f(x_m) + Y_m \ge 0
X_arr{3}=X_consZ; Y_arr{3}=Y_consZ; C_arr{3}=C_consZ; D_arr{3}=D_consZ;%% Equality C^\top f(x_m) + Y_m = 0

% GENERATES GRID AND TRUE VECTOR FIELD
L_grid=linspace(0,1,n_grid)';

stringCommandLeft='['; stringCommandRight2='[';
for i=1:dim_in
    stringCommandLeft=[stringCommandLeft,['x_grid',num2str(i),',']];
    stringCommandRight2=[stringCommandRight2,['x_grid',num2str(i),'(:),']];
end
stringCommandLeft=[stringCommandLeft(1:end-1),']'];
stringCommandRight2=[stringCommandRight2(1:end-1),']'];
stringCommandRight=['= ndgrid(',repmat('L_grid,',1,n_segments),repmat('L_grid,',1,n_segments)]; %2*pi*
stringCommandRight=[stringCommandRight(1:end-1),');'];
eval([stringCommandLeft,stringCommandRight]);
X_grid=eval(stringCommandRight2);
Y_grid=robot_arm(X_grid);

%AUTOMATED PART - FORMATING TO UNIQUE EVALUATION POINTS
[X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq,...
    D_arr_uniq,D_cellArr_uniq_Ineq,bool_arr_C_Ineq_within_C_sol] =...
    ConvertingConsProblem_toUniquePoints_wDeriv(X_arr,Y_arr,C_arr,D_arr,X_grid,Y_grid);
end