function eta_arr = EtaComputingVanLoan(t_arr, delta_arr, nbPointsApprox, A, B, C, ratio)
%ETACOMPUTINGVANLOAN outputs eta_arr (#ineq_constraints x 1), 
% which is a vector of eta parameters 1*#unique_ineq_points. 
% This is achieved by approximating the supremum of the definition of eta_m over the whole interval
% taking the supremum over nbPointsApprox grid samples in the interval [t_arr(m)-delta_arr(m),t_arr(m)+delta_arr(m)]
% The matrices K(s,t_m) are approximated using Van Loan's trick described in "Computing Integrals Involving the
% Matrix Exponential, CHARLES F. VAN LOAN, TAC,1978", then the constraint vectors C are applied.
% In practice very few points nbPointsApprox are required. For refined 
% grids t_arr, two points are often enough, the maximum being often
% attained at the border of [t_arr(m)-delta_arr(m),t_arr(m)+delta_arr(m)].

% Inputs:
% t_arr: 1x#time_points, list of time points where to compute eta_arr
% delta_arr: 1x#time_points, list of "radius" of interval around each time point
% A: NxN, B:NxP, matrices defining x'=Ax+Bu 
% C: cell array of #constraint_points x 1, each cell contains a matrix NxN_cons_m 
% of N_cons_m constraint vectors Nx1 to be applied at point t_m, i.e.
% c_{m,j}^\top x(t_m), j\le N_cons_m.
% ratio: factor in front of the uncontrolled part of the kernel (related to the choice of norm)
% In the experiment the ratio is taken to be 0 since the trajectories start
% at the origin. 

nc=size(C,2); N=size(A,1); nbPointsApprox=ceil(nbPointsApprox/2)*2;
eta_arr=zeros(length(t_arr),nc);
% eta_arr=zeros(length(t_arr),nc,nn);
% delta_arr=delta.*ones(nn,1);
tic
for count=1:length(t_arr)
    nn=nbPointsApprox+1; delta=delta_arr(count);
    listT=t_arr(count)+linspace(-delta,delta,nn);
    BigMat_cellArrKtutu=cell(nn,1); BigMat_cellArrKttu=cell(nn,1);
    temptt = expm([A B*B';zeros(N,N) -A']*t_arr(count));%K(t,t)
    for i=1:nn %Van Loan's technique for computing Gramians
        temp = expm([A B*B';zeros(N,N) -A']*listT(i));
        BigMat_cellArrKtutu{i}=ratio*temp(1:N,1:N)*temp(1:N,1:N)'+temp(1:N,N+1:2*N)*temp(1:N,1:N)';%K(t+u,t+u)
        if i<nn/2
        BigMat_cellArrKttu{i}=ratio*temp(1:N,1:N)*temptt(1:N,1:N)'+temp(1:N,N+1:2*N)*temptt(1:N,1:N)';%K(t+u,t) u<0
        else
        BigMat_cellArrKttu{i}=ratio*temp(1:N,1:N)*temptt(1:N,1:N)'+temp(1:N,1:N)*temptt(1:N,N+1:2*N)'; %K(t+u,t) u>0   
        end
    end
    for i=1:nn
     temp_mat=  BigMat_cellArrKtutu{i} + BigMat_cellArrKtutu{ceil(nn/2)} - 2*BigMat_cellArrKttu{i};
    for j=1:nc
%     eta_arr(count,j,i)= sqrt(abs(C(:,j)'*temp_mat*C(:,j)));
    eta_arr(count,j)= max(eta_arr(count,j),sqrt(abs(C(:,j)'*temp_mat*C(:,j))));
    end
    end
    
end
elapsedTime=toc;
disp(['Van Loan: finished eta ' num2str(elapsedTime) 's']);
end

% nb_exp=40; lambda=1E-2; c_1=[0 -1 0]'; c_2=[0 0 1]'; C=[c_1, c_2]; delta_arr=linspace(1E-5,0.10,nb_exp);
% eta_arr_1 = EtaComputingVanLoan(repmat(0.5,1,nb_exp),delta_arr, 2, A, B, C, lambda);
% eta_arr_2 = EtaComputingVanLoan(repmat(0.9,1,nb_exp), delta_arr, 2, A, B, C, lambda);
% 
% figure
% hold on
% plot(delta_arr, eta_arr_1(:,1), 'b',delta_arr, eta_arr_1(:,2), 'r')
% plot(delta_arr, eta_arr_2(:,1), '--b',delta_arr, eta_arr_2(:,2), '--r')