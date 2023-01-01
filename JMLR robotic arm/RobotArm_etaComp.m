function [eta_mat,D_list,rho_hyp_mat,r_hyp_mat] = RobotArm_etaComp(sig_param,radius,kFuncDeriv_mat,nb_samplesCube,dim_in,D_consIneq)
% ROBOTARM_ETACOMP outputs eta_map (size(kFuncDeriv_mat)), 
% which is a cell array of eta parameters for the various derivatives D of
% the translation-invariant kernel selected when they appear in the constraint of interest. 
% rho_hyp_mat is for the hyperplane covering, r_hyp_mat for ball covering.
% The computations are achieved by approximating the supremum of the definition of eta_m
% taking the supremum over nb_samplesCube samples in the hypersphere. The
% samples are obtained using randsphere.m.

% Inputs:
% sig_param: parameters of the kernel (e.g. the bandwidths)
% radius: the radius of the hyperball over which eta is defined
% kFuncDeriv_mat: the cell array of the matlab functions of the derivatives of
% the kernel
% nb_samplesCube: the number of samples taken to estimate eta
% dim_in: the dimension of the inputs in the kernel
% D_consIneq: the array of indices (0,1,0,..) of the derivatives appearing
% in the constraints.

D_list=unique(full(D_consIneq'),'rows')';
center=zeros(1,dim_in);
eta_mat=cell(size(kFuncDeriv_mat));
rho_hyp_mat=cell(size(kFuncDeriv_mat));
r_hyp_mat=cell(size(kFuncDeriv_mat));

% loc_samples=(rand(nb_samplesCube,dim_in)-ones(1,dim_in)/2)*2*radius+center;
loc_samples=randsphere(nb_samplesCube,dim_in,radius);
for i=1:size(D_list,2)
    temp_idxDeriv=num2cell([D_list(:,i)'+1,D_list(:,i)'+1]);
    temp_DiDik=kFuncDeriv_mat{temp_idxDeriv{:}};
    tmp_valEtai=0;
    tmp_valRhoi=1E6;
for j=1:nb_samplesCube
    tmp_val=temp_DiDik(center,center,sig_param)...
        +temp_DiDik(loc_samples(j,:),loc_samples(j,:),sig_param)...
        -temp_DiDik(center,loc_samples(j,:),sig_param)...
        -temp_DiDik(loc_samples(j,:),center,sig_param);
    tmp_val_hyp=temp_DiDik(loc_samples(j,:),center,sig_param);
    tmp_valEtai=max(tmp_valEtai,sqrt(tmp_val));
    tmp_valRhoi=min(tmp_valRhoi,tmp_val_hyp);
end
    eta_mat{temp_idxDeriv{:}}=tmp_valEtai;
    r_hyp_mat{temp_idxDeriv{:}}=sqrt(temp_DiDik(center,center,sig_param));
    rho_hyp_mat{temp_idxDeriv{:}}=tmp_valRhoi;
end
end