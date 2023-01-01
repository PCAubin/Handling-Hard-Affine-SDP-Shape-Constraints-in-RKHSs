% The code below allow to reproduce a minimal numerical result for the 
% robotic arm experiment of the article
% P.-C. Aubin-Frankowski, Z. Szabo "Handling Hard Affine SDP Shape
% Constraints in RKHSs ", JMLR 2022
% This partial experiment should take about 20s to run on a normal computer.

% The experiment requires having installed YALMIP, we used a free academic
% license for MOSEK, the solver used in this code 
% It can be retrieved at: https://www.mosek.com/downloads/

addpath(genpath('C:\Users\pierr\Documents\MATLAB\YALMIP-master'))
clear all
close all

tic
% CHANGE THE PARAMETERS BELOW TO VARY THE NUMBER OF OBSERVED POINTS AND THE
% NUMBER OF SEGMENTS OF THE ROBOT (EXPERIMENT WILL TAKE ABOUT 20s for
% n_segments=2
% n_segments=1; n_samples=20;
n_segments=2; n_samples=40;
% n_segments=3; n_samples=40;

dim_in=2*n_segments;dim_out=3;%
n_grid=4;n_consIneq_eval=400;n_consIneq_disc=2^dim_in;sigma_noise=0.2;%3^dim_in
load(['GaussianFuncDeriv_mat_dim_',num2str(dim_in),'_Dmax_1.mat']);
% load(['Matern52FuncDeriv_mat_dim_',num2str(dim_in),'_Dmax_1.mat']);

% SET THE BOOLEANS BELOW TO 0 (NO COMPUTATION) OR 1 TO SELECT FOR WHICH
% CONSTRAINT YOU WANT TO BE TESTED
bool_solDisc=1;bool_solSoC_ball=1;bool_solSoC_hyp=0;bool_solSoS=0;
% DELTAFACTOR IS VERY IMPORTANT, IT IS THE INVERSE OF THE RADIUS OF THE BALLS WHERE THE
% SOC CONSTRAINT IS APPLIED, SET IT LARGE TO SATISFY THE NONNEGATIVITY
% CONSTRAINT ONLY VERY LOCALLY, OR SMALL AND YOU RISK GETTING A ZERO
% SOLUTION
deltaFactor=50;

%TO HAVE A GRID FOR EVALUATING THE CONSTRAINT VIOLATION CHANGE BELOW,
%OTHERWISE POINTS ARE SAMPLED RANDOMLY IN THE HYPERCUBE. HIGH RISK OF
%STRANGE RESULTS (NULL VIOLATION) DUE TO PERIODICITY IF GRID IS APPLIED
boolGrid_consIneq_disc=0;
if boolGrid_consIneq_disc
    n_consIneq_disc=200;
end
% ONLY WORKS FOR ONE SINGLE PARTIAL DERIVATIVE IN EACH CONSTRAINT/NOT GENERAL DIFF OPERATORS
[X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq,...
    D_arr_uniq,D_cellArr_uniq_Ineq,X_consIneq_eval,C_consIneq_eval,D_consIneq_eval,bool_arr_C_Ineq_within_C_sol] =...
    RobotArm_generation(n_grid,n_samples,n_consIneq_eval,n_consIneq_disc,...%
    sigma_noise,n_segments,boolGrid_consIneq_disc);%_oneSidedCons
elapsedTime=toc;
disp(['Generating points time:' num2str(elapsedTime) 's'])

%%
%Agrell's params: [0.60671424 0.80578935 0.14750236 0.16083182] 
sig_param=[ones(1,n_segments), 0.2*ones(1,n_segments)]; 
sigSoS=0.2;%sigK/sqrt(2); 
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
% Sig_corelMat=eye(dim_out);

[mat_obj,mat_consIneq,mat_consEq,mat_consIneq_eval,mat_Grid,mat_norm,R_SoS,lambdaK,...
    v_field_naive,valsIneq_field_naive_eval]...
     = RobotArm_EvalKernel(X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq, D_arr_uniq,...
    kFuncDeriv_mat,sig_param,Sig_corelMat,lambdaK,hSoS,X_consIneq_eval,C_consIneq_eval,D_consIneq_eval);
elapsedTime=toc;
disp(['Computing kernel matrices time:' num2str(elapsedTime) 's'])
%

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

disp(['Average of eta: '...
    num2str(mean(etaK_arr)) '. Average of rho: ' num2str(mean(rhoK_arr)) '.'])

blobCheckGeom=etaKSig_arr.^2-2*(rKSig_arr.^2-rhoKSig_arr);
blobRatio=CSig_arr.*rhoK_arr./rK_arr;
%

etaFactor=1;%A_solSoC_hyp
rhoFactor_hyp=1;%0.99/max(rhoK_arr);
[A_solDisc,A_solSoC_ball,A_solSoC_hyp,A_solSoS] = RobotArm_solver_YALMIP(...
    Y_arr_uniq,mat_obj,mat_consIneq,...
    mat_consEq,mat_norm,R_SoS,lambdaK,lambdaSoS,etaKSig_arr*etaFactor,CSig_arr.*rhoK_arr*rhoFactor_hyp,rK_arr,...%CSig_arr.*%rK_arr*rhoFactor_hyp
    bool_solDisc,bool_solSoC_ball,bool_solSoC_hyp,bool_solSoS,bool_arr_C_Ineq_within_C_sol);
%
if ~bool_solSoC_hyp
    A_solSoC_hyp=A_solSoC_ball;
end
if ~bool_solSoS
    A_solSoS=zeros(size(A_solSoC_ball));
end
% A_solSoS=zeros(size(A_solDisc));
v_field_Ineq=reshape((mat_Grid*A_solDisc)',dim_out,[])';
v_field_Ineq_SoS=reshape((mat_Grid*A_solSoS)',dim_out,[])';
v_field_Ineq_SoC_ball=reshape((mat_Grid*A_solSoC_ball)',dim_out,[])';
v_field_Ineq_SoC_hyp=reshape((mat_Grid*A_solSoC_hyp)',dim_out,[])';

valsIneq_field_Ineq_eval=mat_consIneq_eval*A_solDisc;
valsIneq_field_Ineq_SoS_eval=mat_consIneq_eval*A_solSoS;
valsIneq_field_Ineq_SoC_ball_eval=mat_consIneq_eval*A_solSoC_ball;
valsIneq_field_Ineq_SoC_hyp_eval=mat_consIneq_eval*A_solSoC_hyp;

v_field_true=Y_grid;
error_arr_naive = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_naive,valsIneq_field_naive_eval);
error_arr_Ineq = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_Ineq,valsIneq_field_Ineq_eval);
error_arr_Ineq_SoS = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_Ineq_SoS,valsIneq_field_Ineq_SoS_eval);
error_arr_Ineq_SoC_ball = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_Ineq_SoC_ball,valsIneq_field_Ineq_SoC_ball_eval);
error_arr_Ineq_SoC_hyp = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,v_field_Ineq_SoC_hyp,valsIneq_field_Ineq_SoC_hyp_eval);
error_arr_Ineq_fail = RobotArm_ErrorComputationVectorField_wBdary(v_field_true,zeros(size(v_field_naive)),zeros(size(valsIneq_field_naive_eval)));
error_mat=[error_arr_naive;error_arr_Ineq;error_arr_Ineq_SoC_ball;error_arr_Ineq_SoC_hyp;error_arr_Ineq_SoS;error_arr_Ineq_fail]%
%THIS CORRESPONDS TO THE ERROR [Q2,L2Err,LinftyErr,BdryErrLinfty,BdryErrL1]
%FOR THE COLUMNS, AND FOR THE ROWS: NO-CONS, DISC CONS, BALL SOC CONS, HYP+BALL SOC CONS,
%KSOS CONS, AND THE FAILURE CASE (SOLUTION IS NULL, HENCE TRIVIAL)

elapsedTime=toc;
disp(['Total time:' num2str(elapsedTime) 's'])

% %% VARIOUS PLOTTING FUNCTION
% figure
% hold on
% % scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq(:,1),'b')   
% % scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoS(:,1),'m')  
% % scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoC_ball(:,1),'g')  
% scatter3(X_consIneq(:,1),X_consIneq(:,2), zeros(size(X_consIneq,1),1),'*k')  
% scatter3(X_consIneq_eval(:,1),X_consIneq_eval(:,2), valsIneq_field_Ineq_SoC_ball_eval,'og')
% scatter3(X_consIneq_eval(:,1),X_consIneq_eval(:,2), valsIneq_field_Ineq_SoC_hyp_eval,'or')
% scatter3(X_consIneq_eval(:,1),X_consIneq_eval(:,2), valsIneq_field_Ineq_eval,'ob')
% % scatter3(X_samples(:,1),X_samples(:,2), Y_samples_vec(:,1),'r')
% % scatter3(X_grid(:,1),X_grid(:,2), v_field_true(:,1),'k')  
% % scatter3(X_grid(:,1),X_grid(:,2), v_field_naive(:,1),'c')  
% % axis tight 
% % axis equal
% view(90,0)
% title('Plot constraints')
% 
% %% PLOTS
% % v_field_Ineq=reshape((mat_Grid*A_solDisc)',dim_out,[])';
% % v_field_true=Y_grid;
% %isequal(Y_samples_vec,robot_arm(X_samples))
% figure
% hold on
% scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq(:,1),'b')   
% % scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoS(:,1),'m')  
% scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoC_ball(:,1),'g')  
% scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoC_hyp(:,1),'r')  
% scatter3(X_consIneq(:,1),X_consIneq(:,2), zeros(size(X_consIneq,1),1),'*k')  
% % scatter3(X_consIneq_eval(:,1),X_consIneq_eval(:,2), zeros(size(X_consIneq_eval,1),1),'ok')
% % scatter3(X_samples(:,1),X_samples(:,2), Y_samples_vec(:,1),'r')
% scatter3(X_grid(:,1),X_grid(:,2), v_field_true(:,1),'k')  
% % scatter3(X_grid(:,1),X_grid(:,2), v_field_naive(:,1),'c')  
% % axis tight 
% % axis equal
% view(90,0)
% title('1st component')
% %%
% figure
% hold on
% scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq(:,2),'b')   
% % scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoS(:,2),'m')  
% scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoC_ball(:,2),'g')  
% scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoC_hyp(:,2),'r')  
% % scatter3(X_consIneq(:,1),X_consIneq(:,2), zeros(size(X_consIneq,1),1),'*k')  
% % scatter3(X_consIneq_eval(:,1),X_consIneq_eval(:,2), zeros(size(X_consIneq_eval,1),1),'ok')
% % scatter3(X_samples(:,1),X_samples(:,2), Y_samples_vec(:,2),'r')
% scatter3(X_grid(:,1),X_grid(:,2), v_field_true(:,2),'k')  
% % axis tight 
% % axis equal
% view(90,0)
% title('2nd component')
% %% 
% figure
% hold on
% scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq(:,3),'b')   
% % scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoS(:,3),'m')  
% scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoC_ball(:,3),'g')  
% scatter3(X_grid(:,1),X_grid(:,2), v_field_Ineq_SoC_hyp(:,3),'r')  
% % scatter3(X_consIneq_eval(:,1),X_consIneq_eval(:,2), zeros(size(X_consIneq_eval,1),1),'ok')
% % scatter3(X_samples(:,1),X_samples(:,2), Y_samples_vec(:,3),'r')
% scatter3(X_grid(:,1),X_grid(:,2), v_field_true(:,3),'k')  
% % axis tight 
% % axis equal
% view(90,0)
% title('3rd component')