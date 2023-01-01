% The code below allow to reproduce a table numerical result for the 
% robotic arm experiment of the article
% P.-C. Aubin-Frankowski, Z. Szabo "Handling Hard Affine SDP Shape
% Constraints in RKHSs ", JMLR 2022
% This partial experiment should take about 40s to run on a normal computer.

% The experiment requires having installed YALMIP, we used a free academic
% license for MOSEK, the solver used in this code 
% It can be retrieved at: https://www.mosek.com/downloads/

addpath(genpath('C:\Users\pierr\Documents\MATLAB\YALMIP-master'))
clear all
close all

% n_segments=1; n_samples=20;
n_segments=2; n_samples=40;
% n_segments=3; n_samples=40;
dim_in=2*n_segments;dim_out=3;%

%PLACE TO CHANGE THE PARAMETERS
deltaFactor_arr=[50 150 200];%[5 10:2:20, 30:10:50, 70, 100];
for nb_repeats=1%:20%18
for nbDisc_pts_perSide=3%:5
startTimeGen=tic;

n_grid=4;n_consIneq_eval=500;n_consIneq_disc=nbDisc_pts_perSide^dim_in;sigma_noise=0.2;%3^dim_in
load(['GaussianFuncDeriv_mat_dim_',num2str(dim_in),'_Dmax_1.mat']);
% load(['Matern52FuncDeriv_mat_dim_',num2str(dim_in),'_Dmax_1.mat']);

boolGrid_consIneq_disc=0;
% ONLY WORKS FOR ONE SINGLE PARTIAL DERIVATIVE IN EACH CONSTRAINT/NOT GENERAL DIFF OPERATORS
[X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq,...
    D_arr_uniq,D_cellArr_uniq_Ineq,X_consIneq_eval,C_consIneq_eval,D_consIneq_eval,bool_arr_C_Ineq_within_C_sol] =...
    RobotArm_generation(n_grid,n_samples,n_consIneq_eval,n_consIneq_disc,...%
    sigma_noise,n_segments,boolGrid_consIneq_disc);%_oneSidedCons
elapsedTimeGen=toc(startTimeGen);
disp(['Generating points time:' num2str(elapsedTimeGen) 's'])
%%
startTimeKer=tic;
sig_param=[ones(1,n_segments), 0.2*ones(1,n_segments)]; sigSoS=0.2;%sig_param(end)/sqrt(2); 
hgauss = @(u,sig) exp(-u.^2/(2*sig^2)); hcauchy = @(u,sig) 1./(1+u.^2/sig^2);
hSoS=@(u) hcauchy(u,sigSoS);
lambdaK=1E-6/sqrt(n_samples); lambdaSoS=1E-8;
% lambdaK=1E-8/sqrt(n_samples); % if dim=2

X_samples=X_arr_uniq{2}; Y_samples=Y_arr_uniq{2}; Y_samples_vec=full(reshape(sum(Y_samples),dim_out,[]))';
X_grid=X_arr_uniq{end}; Y_grid=Y_arr_uniq{end};
C_sol=C_arr_uniq{1}; 
X_consIneq=X_arr_uniq{3}; C_consIneq=C_arr_uniq{3}; D_consIneq=D_arr_uniq{3};

Sig_corelMat_th=RobotArm_estimateCovOutput(dim_in);
Sig_corelMat_emp=cov(Y_samples_vec);
Sig_corelMat=Sig_corelMat_th;

[mat_obj,mat_consIneq,mat_consEq,mat_consIneq_eval,mat_Grid,mat_norm,R_SoS,lambdaK,...
    v_field_naive,valsIneq_field_naive_eval]...
     = RobotArm_EvalKernel(X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq, D_arr_uniq,...
    kFuncDeriv_mat,sig_param,Sig_corelMat,lambdaK,hSoS,X_consIneq_eval,C_consIneq_eval,D_consIneq_eval);
elapsedTimeKer=toc(startTimeKer);
% disp(['Computing kernel matrices time:' num2str(elapsedTimeKer) 's'])


for deltaFactor= deltaFactor_arr
bool_solSoC_ball=1; bool_solDisc=1;bool_solSoC_hyp=1; bool_solSoS=1;
if nbDisc_pts_perSide > 1
	bool_solSoS=0;
end

if nbDisc_pts_perSide > 2
	bool_solSoC_hyp=0;
end

if bool_solSoC_ball||bool_solSoC_hyp
nbPts_per_dimCons=max(2,floor(nthroot(n_consIneq_disc,dim_in)));
    deltaDesi=1/nbPts_per_dimCons/2;
    deltaDesi=deltaDesi/deltaFactor;
[eta_mat,D_list,rho_hyp_mat,r_hyp_mat] = RobotArm_etaComp(sig_param,deltaDesi,kFuncDeriv_mat,1000,dim_in,D_consIneq);
    temp_cellSig_corelMat = repmat({Sig_corelMat},1,size(C_sol,1)/size(Sig_corelMat,1));
    CSig_arr=sqrt(diag(C_consIneq'*blkdiag(temp_cellSig_corelMat{:})*C_consIneq));
    etaK_arr=zeros(size(CSig_arr)); rhoK_arr=zeros(size(CSig_arr)); rK_arr=zeros(size(CSig_arr));
    for i=1:size(C_consIneq,2)
    temp_idxDeriv=num2cell([D_consIneq(:,i)'+1,D_consIneq(:,i)'+1]);
        etaK_arr(i)=eta_mat{temp_idxDeriv{:}};
        rhoK_arr(i)=rho_hyp_mat{temp_idxDeriv{:}};
        rK_arr(i)=r_hyp_mat{temp_idxDeriv{:}};
    end%rho_hyp_mat,r_hyp_mat
%     idx_neg = any(C_consIneq <0,1);
%     etaK_arr(idx_neg,:) = 0;
%     rhoK_arr(idx_neg,:) = 0;
%     rK_arr(idx_neg,:) = 0;
    etaKSig_arr=CSig_arr.*etaK_arr;
    rhoKSig_arr=CSig_arr.^2.*rhoK_arr;
    rKSig_arr=CSig_arr.*rK_arr;
else
eta_mat=NaN;
end

% disp(['Average of eta: '...
%     num2str(mean(etaK_arr)) '. Average of rho: ' num2str(mean(rhoK_arr)) '.'])

etaFactor=1;%A_solSoC_hyp
rhoFactor_hyp=1;%0.99/max(rhoK_arr);
[A_solDisc,A_solSoC_ball,A_solSoC_hyp,A_solSoS,time_matrix] = RobotArm_solver_YALMIP(...
    Y_arr_uniq,mat_obj,mat_consIneq,...
    mat_consEq,mat_norm,R_SoS,lambdaK,lambdaSoS,etaKSig_arr*etaFactor,CSig_arr.*rhoK_arr*rhoFactor_hyp,rK_arr,...%CSig_arr.*%rK_arr*rhoFactor_hyp
    bool_solDisc,bool_solSoC_ball,bool_solSoC_hyp,bool_solSoS,bool_arr_C_Ineq_within_C_sol);

v_field_Ineq=reshape((mat_Grid*A_solDisc)',dim_out,[])';
valsIneq_field_Ineq_eval=mat_consIneq_eval*A_solDisc;
v_field_true=Y_grid;
error_arr_naive = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_naive,valsIneq_field_naive_eval);
error_arr_Ineq = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_Ineq,valsIneq_field_Ineq_eval);
error_arr_Ineq_fail = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,zeros(size(v_field_naive)),zeros(size(valsIneq_field_naive_eval)));
% error_mat=[error_arr_naive;error_arr_Ineq;error_arr_Ineq_SoC_ball;error_arr_Ineq_SoC_hyp;error_arr_Ineq_SoS;error_arr_Ineq_fail]%

elapsedTime=toc(startTimeGen);
disp(['Total time of experiment #' num2str(nb_repeats) ':' num2str(elapsedTime) 's for '...
    num2str(nbDisc_pts_perSide) 'samples per side and deltaFactor=' num2str(deltaFactor)])

saveFileName=['table_results_RobotArm_' num2str(dim_in) 'D_fullGaussian_wErrBdary_' num2str(n_consIneq_eval) 'ptEval_dummy.txt'];
fileID = fopen(saveFileName,'a');

fprintf(fileID,'%d \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t%d \t %d \t %d \t %.3e \t %.3e \t %.3e  \t %.3e  \t %.3e \t %.1e \t %.1e \t %.1e\t %.1e\n',...
        0, error_arr_naive, nbDisc_pts_perSide, size(D_consIneq,2), deltaFactor, sigma_noise, sigSoS, sig_param(end), lambdaK, lambdaSoS, time_matrix(1,:),elapsedTime);
    
fprintf(fileID,'%d \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t%d \t %d \t %d \t %.3e \t %.3e \t %.3e  \t %.3e  \t %.3e \t %.1e \t %.1e \t %.1e\t %.1e\n',...
        1, error_arr_Ineq, nbDisc_pts_perSide,size(D_consIneq,2), deltaFactor, sigma_noise, sigSoS, sig_param(end), lambdaK, lambdaSoS, time_matrix(2,:),elapsedTime);
if bool_solSoC_ball   
    valsIneq_field_Ineq_SoC_ball_eval=mat_consIneq_eval*A_solSoC_ball;
    v_field_Ineq_SoC_ball=reshape((mat_Grid*A_solSoC_ball)',dim_out,[])';
    error_arr_Ineq_SoC_ball = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_Ineq_SoC_ball,valsIneq_field_Ineq_SoC_ball_eval);
    fprintf(fileID,'%d \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t%d \t %d \t %d \t %.3e \t %.3e \t %.3e  \t %.3e  \t %.3e \t %.1e \t %.1e \t %.1e\t %.1e\n',...
        2, error_arr_Ineq_SoC_ball, nbDisc_pts_perSide,size(D_consIneq,2), deltaFactor, sigma_noise, sigSoS, sig_param(end), lambdaK, lambdaSoS, time_matrix(3,:),elapsedTime);
end    
if bool_solSoC_hyp
    valsIneq_field_Ineq_SoC_hyp_eval=mat_consIneq_eval*A_solSoC_hyp;    
    v_field_Ineq_SoC_hyp=reshape((mat_Grid*A_solSoC_hyp)',dim_out,[])';
    error_arr_Ineq_SoC_hyp = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_Ineq_SoC_hyp,valsIneq_field_Ineq_SoC_hyp_eval);
fprintf(fileID,'%d \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t%d \t %d \t %d \t %.3e \t %.3e \t %.3e  \t %.3e  \t %.3e \t %.1e \t %.1e \t %.1e\t %.1e\n',...
        3, error_arr_Ineq_SoC_hyp, nbDisc_pts_perSide,size(D_consIneq,2), deltaFactor, sigma_noise, sigSoS, sig_param(end), lambdaK, lambdaSoS, time_matrix(4,:),elapsedTime);
end
if bool_solSoS
    valsIneq_field_Ineq_SoS_eval=mat_consIneq_eval*A_solSoS;
    v_field_Ineq_SoS=reshape((mat_Grid*A_solSoS)',dim_out,[])';
    error_arr_Ineq_SoS = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_Ineq_SoS,valsIneq_field_Ineq_SoS_eval);
fprintf(fileID,'%d \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t%d \t %d \t %d \t %.3e \t %.3e \t %.3e  \t %.3e  \t %.3e \t %.1e \t %.1e \t %.1e\t %.1e\n',...
        4, error_arr_Ineq_SoS, nbDisc_pts_perSide,size(D_consIneq,2), deltaFactor, sigma_noise, sigSoS, sig_param(end), lambdaK, lambdaSoS, time_matrix(5,:),elapsedTime);
end
fprintf(fileID,'%d \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t%d \t %d \t %d \t %.3e \t %.3e \t %.3e  \t %.3e  \t %.3e \t %.1e \t %.1e \t %.1e\t %.1e\n',...
        5, error_arr_Ineq_fail, nbDisc_pts_perSide,size(D_consIneq,2), deltaFactor, sigma_noise, sigSoS, sig_param(end), lambdaK, lambdaSoS, time_matrix(6,:),elapsedTime);
fclose(fileID); 
end
end
end