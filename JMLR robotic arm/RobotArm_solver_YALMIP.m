function [A_solDisc,A_solSoC_ball,A_solSoC_hyp,A_solSoS,time_matrix] = RobotArm_solver_YALMIP_alternative(...
    Y_arr_uniq,mat_obj,mat_consIneq,...
    mat_consEq,mat_norm,R_SoS,lambdaK,lambdaSoS,etaKSig_arr,rhoKSig_arr,rKSig_arr,...
    bool_solDisc,bool_solSoC_ball,bool_solSoC_hyp,bool_solSoS,bool_arr_C_Ineq_within_C_sol)
%   RobotArm_SOLVER Summary of this function goes here
%   Detailed explanation goes here

% blob=etaKSig_arr.^2-2*(rKSig_arr.^2-rhoKSig_arr);
Y_samples=Y_arr_uniq{2}; Y_consIneq=Y_arr_uniq{3}; Y_consEq=Y_arr_uniq{4};
% dim_in=size(X_sol,2); dim_out=size(Y_grid,2);
% Y_samples_vec=full(reshape(sum(Y_samples),dim_out,[]))';

time_matrix=zeros(6,3);

nb_coeff=size(mat_obj,2); nb_samples=size(Y_samples,1); n_SoS=size(R_SoS,1);
tol_SOC=1E-15;
GnormSol=chol((mat_norm+mat_norm')/2+tol_SOC*eye(size(mat_norm,1)));

A_solDisc=NaN; A_solSoC_ball=NaN; A_solSoC_hyp=NaN; A_solSoS=NaN; 
ops = sdpsettings('solver','mosek','verbose',0,'debug',0);
yalmip('clear')
tic
A = sdpvar(nb_coeff,1);
sdpvar dummyNorm dummyObj

if bool_solDisc
    tic
F = [cone(full(sum(Y_samples))'-mat_obj*A,dummyObj),...
    cone(sqrt(lambdaK*nb_samples)*GnormSol*A,dummyNorm),...
    0 <= full(sum(Y_consIneq))'+ mat_consIneq*A,...
    0 == full(sum(Y_consEq))'+ mat_consEq*A];
optimize(F,dummyObj + dummyNorm,ops);    
A_solDisc=value(A);
elapsedTime=toc;
time_matrix(2,:)=[ans.yalmiptime ans.solvertime elapsedTime];
end

yalmip('clear')

A = sdpvar(nb_coeff,1);
sdpvar dummyNorm dummyObj

if bool_solSoC_ball
tic
F = [cone(full(sum(Y_samples))'-mat_obj*A,dummyObj),...
    cone(sqrt(lambdaK*nb_samples)*GnormSol*A,dummyNorm),...
    etaKSig_arr*norm(GnormSol*A,2) <= full(sum(Y_consIneq))'+ mat_consIneq*A,...
    0 == full(sum(Y_consEq))'+ mat_consEq*A];
optimize(F,dummyObj + dummyNorm,ops);    
A_solSoC_ball=value(A);
elapsedTime=toc;
time_matrix(3,:)=[ans.yalmiptime ans.solvertime elapsedTime];
end

yalmip('clear')

A = sdpvar(nb_coeff,1);
sdpvar dummyNorm dummyObj

if bool_solSoC_hyp
tic
    nb_hyper=size(mat_consIneq,1);
    boolMat_hyp=diag(bool_arr_C_Ineq_within_C_sol);
    boolMat_hyp(:,all(boolMat_hyp == 0))=[];
    boolMat_hyp=double(boolMat_hyp);
B = sdpvar(nb_hyper,1);
F = [cone(full(sum(Y_samples))'-mat_obj*A,dummyObj),...
    cone(sqrt(lambdaK*nb_samples)*GnormSol*A,dummyNorm),...
    0 == full(sum(Y_consEq))'+ mat_consEq*A];
% for i=1:nb_hyper
%     if rKSig_arr(i)~=0
%      F = [F,cone(GnormSol*(A-boolMat_hyp(:,i)*B(i)),B(i)*rhoKSig_arr(i)/rKSig_arr(i))];%rhoKSig_arr(i);
%     end
% end 
F = [F,cone([(B.*rhoKSig_arr./rKSig_arr)';GnormSol*(repmat(A,1,nb_hyper)-boolMat_hyp.*(boolMat_hyp*repmat(B,1,nb_hyper)))])];%rhoKSig_arr(i);
% 
optimize(F,dummyObj + dummyNorm,ops);    
A_solSoC_hyp=value(A);
% disp(['min value of hyp-covering on discretization: ' num2str(min(mat_consIneq*A_solSoC_hyp))])
elapsedTime=toc;
time_matrix(4,:)=[ans.yalmiptime ans.solvertime elapsedTime];
end

yalmip('clear')

A = sdpvar(nb_coeff,1);
sdpvar dummyNorm dummyObj

ops = sdpsettings('solver','mosek','verbose',1,'debug',1);
if bool_solSoS
    tic
    eq_vect=full(sum(Y_consIneq))';
    B = sdpvar(n_SoS,n_SoS);
    F = [cone(full(sum(Y_samples))'-mat_obj*A,dummyObj),...
    cone(sqrt(lambdaK*nb_samples)*GnormSol*A,dummyNorm),...
    0 == full(sum(Y_consEq))'+ mat_consEq*A, B >= 0];
    for i=1:n_SoS
        cons_temp=R_SoS(:,i)*R_SoS(:,i)';
        F = [F,trace(B*cons_temp) == eq_vect(i) + mat_consIneq(i,:)*A];
    end
%     ops = sdpsettings('solver','mosek','verbose',1,'debug',1);
optimize(F,dummyObj + dummyNorm+sqrt(lambdaSoS)*trace(B),ops);    
A_solSoS=value(A);

elapsedTime=toc;
time_matrix(5,:)=[ans.yalmiptime ans.solvertime elapsedTime];
end 


end