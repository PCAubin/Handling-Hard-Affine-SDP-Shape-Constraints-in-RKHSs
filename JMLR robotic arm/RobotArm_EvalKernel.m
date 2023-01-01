function [mat_obj,mat_consIneq,mat_consEq,mat_consIneq_eval,mat_Grid,mat_norm,R_SoS,lambdaK,...
    v_field_naive,valsIneq_field_naive_eval]...
     = RobotArm_EvalKernel(X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq, D_arr_uniq,...
    kFuncDeriv_mat,sig_param,Sig_corelMat,lambdaK,hSoS,X_consIneq_eval,C_consIneq_eval,D_consIneq_eval)
%ROBOTARM_EVALKERNEL computes the various Gram matrices which appear when
%expressing the finite-dimensional solution of the optimization problem
%with discretized (SOC or not) constraints. 
%mat_obj is the Gram matrix of the objective, mat_consIneq of the
%inequality constraints at the sampled points and mat_consEq of the
%equality constraints. mat_consIneq_eval is for the points where the
%inequality constraints are evaluated (L1 violation), mat_Grid for the grid points where
%the function is evaluated (L2 reconstruction error), mat_Grid is the for
%the square norm regularization term, R_SoS is the Cholesky decomposition
%of the kernel SoS Gram matrix, lambdaK is the regularization parameter of
%the unconstrained problem.
% Inputs:
% (X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq, D_arr_uniq: cell array of block diagonal 
% matrices obtained as output of RobotArm_generation
% kFuncDeriv_mat: derivatives in a matlab function form of the base kernel
% k with parameters sig_param
% Sig_corelMat: correlation matrix used to define the matrix-valued kernel
% hSoS: tranlation invariant kernel used for kSoS
% X_consIneq_eval,C_consIneq_eval,D_consIneq_eval: copies of the data, of
% the constraint vector and derivative index required for computing the
% Gram matrices

X_grid=X_arr_uniq{end}; Y_grid=Y_arr_uniq{end};
X_sol=X_arr_uniq{1}; C_sol=C_arr_uniq{1}; D_sol=D_arr_uniq{1};
X_samples=X_arr_uniq{2}; Y_samples=Y_arr_uniq{2}; C_samples=C_arr_uniq{2}; D_samples=D_arr_uniq{2};% objective
X_consIneq=X_arr_uniq{3}; Y_consIneq=Y_arr_uniq{3}; C_consIneq=C_arr_uniq{3}; D_consIneq=D_arr_uniq{3};% inequality
X_consEq=X_arr_uniq{4}; Y_consEq=Y_arr_uniq{4}; C_consEq=C_arr_uniq{4}; D_consEq=D_arr_uniq{4};% equality

dim_in=size(X_sol,2); dim_out=size(Y_grid,2);
Y_samples_vec=full(reshape(sum(Y_samples),dim_out,[]))';

%Build SoS matrix
tol=1E-8;
vec_repeat=cell2mat(cellfun(@size,C_cellArr_uniq_Ineq,'UniformOutput',false));
vec_repeat=vec_repeat(:,2);
X_SoS=repelem(X_consIneq,vec_repeat,1);
n_SoS=size(X_SoS,1);
GX_SoS=hSoS(squareform(pdist(X_SoS,'euclidean')));
R_SoS=chol(GX_SoS+tol*eye(n_SoS));

GXX_C=zeros(size(C_sol,2));
GXsamplesX_C=zeros(size(C_samples,2),size(C_sol,2));
GXIneqX_C=zeros(size(C_consIneq,2),size(C_sol,2));
GXEqX_C=zeros(size(C_consEq,2),size(C_sol,2));

GXgridX_C=zeros(size(X_grid,1),size(C_sol,2));
GXevalX_C=zeros(size(X_consIneq_eval,1),size(C_sol,2));

for j=1:size(C_sol,2)
    temp_idx_xj=floor((find(C_sol(:,j),1)-1)/dim_out)+1;
    for i=1:size(C_sol,2) 
    temp_idx_xi=floor((find(C_sol(:,i),1)-1)/dim_out)+1;
    temp_idxDeriv=num2cell([D_sol(:,i)'+1,D_sol(:,j)'+1]);
    GXX_C(i,j)=kFuncDeriv_mat{temp_idxDeriv{:}}(X_sol(temp_idx_xi,:),X_sol(temp_idx_xj,:),sig_param);
    end
    for i=1:size(C_samples,2) 
    temp_idx_xi=floor((find(C_samples(:,i),1)-1)/dim_out)+1;
    temp_idxDeriv=num2cell([D_samples(:,i)'+1,D_sol(:,j)'+1]);
    GXsamplesX_C(i,j)=kFuncDeriv_mat{temp_idxDeriv{:}}(X_sol(temp_idx_xi,:),X_sol(temp_idx_xj,:),sig_param);
    end
    for i=1:size(C_consIneq,2) 
    temp_idx_xi=floor((find(C_consIneq(:,i),1)-1)/dim_out)+1;
    temp_idxDeriv=num2cell([D_consIneq(:,i)'+1,D_sol(:,j)'+1]);
    GXIneqX_C(i,j)=kFuncDeriv_mat{temp_idxDeriv{:}}(X_sol(temp_idx_xi,:),X_sol(temp_idx_xj,:),sig_param);
    end
    for i=1:size(C_consEq,2) 
    temp_idx_xi=floor((find(C_consEq(:,i),1)-1)/dim_out)+1;
    temp_idxDeriv=num2cell([D_consEq(:,i)'+1,D_sol(:,j)'+1]);
    GXEqX_C(i,j)=kFuncDeriv_mat{temp_idxDeriv{:}}(X_sol(temp_idx_xi,:),X_sol(temp_idx_xj,:),sig_param);
    end
    for i=1:size(X_grid,1)   
    temp_idxDeriv=num2cell([zeros(1,dim_in)+1,D_sol(:,j)'+1]);
    GXgridX_C(i,j)=kFuncDeriv_mat{temp_idxDeriv{:}}(X_grid(i,:),X_sol(temp_idx_xj,:),sig_param);
    end
    for i=1:size(X_consIneq_eval,1)   
    temp_idxDeriv=num2cell([D_consIneq_eval(i,:)+1,D_sol(:,j)'+1]);
    GXevalX_C(i,j)=kFuncDeriv_mat{temp_idxDeriv{:}}(X_consIneq_eval(i,:),X_sol(temp_idx_xj,:),sig_param);
    end
end
%EXCLUSIVELY FOR MATRIX-VALUED KERNELS THAT ARE DECOMPOSABLE (of the form k(x,y)*Sigma)
BigSig_corelMat=repmat(Sig_corelMat,size(X_sol,1),size(X_sol,1));

temp_matProd= BigSig_corelMat*C_sol;

mat_Grid=repelem(GXgridX_C,repmat(dim_out,1,size(GXgridX_C,1)),1).*...
    (repmat(Sig_corelMat,size(X_grid,1),size(X_sol,1))*C_sol);
mat_obj= GXsamplesX_C.*(C_samples'*temp_matProd);%C_samples'*temp_matProd*C_sol; %
mat_consIneq=GXIneqX_C.*(C_consIneq'*temp_matProd);
mat_consEq=GXEqX_C.*(C_consEq'*temp_matProd);
mat_norm=GXX_C.*(C_sol'*temp_matProd);
mat_consIneq_eval=GXevalX_C.*(C_consIneq_eval*repmat(Sig_corelMat,1,size(X_sol,1))*C_sol);
%FI

%CAREFUL: HARD CODED FOR ZERO DERIVATIVES IN THE SAMPLES
G_kernel1D_XsamplesXsamples=zeros(size(X_samples,1));
G_kernel1D_XgridXsamples=zeros(size(X_grid,1),size(X_samples,1));
G_kernel1D_XevalXsamples=zeros(size(X_consIneq_eval,1),size(X_samples,1));

for j=1:size(X_samples,1)
    for i=1:size(X_samples,1)   
    temp_idxDeriv=num2cell([zeros(1,dim_in)+1,zeros(1,dim_in)+1]);
    G_kernel1D_XsamplesXsamples(i,j)=kFuncDeriv_mat{temp_idxDeriv{:}}(X_samples(i,:),X_samples(j,:),sig_param);
    end
    for i=1:size(X_grid,1)   
    temp_idxDeriv=num2cell([zeros(1,dim_in)+1,zeros(1,dim_in)+1]);
    G_kernel1D_XgridXsamples(i,j)=kFuncDeriv_mat{temp_idxDeriv{:}}(X_grid(i,:),X_samples(j,:),sig_param);
    end
    for i=1:size(X_consIneq_eval,1)   
    temp_idxDeriv=num2cell([D_consIneq_eval(i,:)+1,zeros(1,dim_in)+1]);
    G_kernel1D_XevalXsamples(i,j)=kFuncDeriv_mat{temp_idxDeriv{:}}(X_consIneq_eval(i,:),X_samples(j,:),sig_param);
    end
end

GXsamplesXsamples=kron(G_kernel1D_XsamplesXsamples,Sig_corelMat);
GXgridXsamples=kron(G_kernel1D_XgridXsamples,Sig_corelMat);
GXevalXsamples=kron(G_kernel1D_XevalXsamples,Sig_corelMat);

% ATTEMPT AT USING gcv.m TO SET AUTOMATICALLY SET lambdaK WHICH IS THE
% REGULARIZATION PARAMETER IN THE UNCONSTRAINED CASE
% [U,S,~] = svd(GXsamplesXsamples);
% [reg_min,~,~] = gcv(U,diag(S),reshape(Y_samples_vec'+1E-1*rand(size(Y_samples_vec')),[],1));%

%    minimize(sum_square(full(sum(Y_samples))'-mat_obj*A)/nb_samples+lambdaK*quad_form(C_sol*A,GXX)) % SQUARE FORM OBJECTIVE

% reg_min=fminbnd(@(l)mygcv1(GXsamplesXsamples,full(sum(Y_samples))',l),1E-15,1E-3); 
% if reg_min>1E-4
%     lambdaK=reg_min;
% else
%     disp("gcv failed")
% end
alpha_naive=(GXsamplesXsamples+...
    eye(size(full(sum(Y_samples))',1))*lambdaK*size(X_samples,1))\(full(sum(Y_samples))');
v_field_naive = reshape((GXgridXsamples*alpha_naive)',dim_out,[])';
valsIneq_field_naive_eval = sum(C_consIneq_eval.*reshape((GXevalXsamples*alpha_naive)',dim_out,[])',2);

end



