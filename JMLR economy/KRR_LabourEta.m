function [etaDeriv10,etaDeriv01,etaDerivConv] = KRR_LabourEta(sigma,deltaDesi,kFuncDeriv_mat,nb_samples,nb_angles,bool_conv)
% KRR_LabourEta outputs three (float) eta values, etaDeriv10,etaDeriv01 for the
% derivatives in x1 and x2, etaDerivConv for the second-order convexity
% constraint.
% The computations are achieved by approximating the supremum in the definition of eta
% taking the supremum over nb_samples randoms samples in square, and nb_angles 
% deterministic angles for the SDP constraints.

% Inputs:
% sigma: parameter of the kernel (e.g. the bandwidth)
% deltaDesi: half the side of the cube over which eta is defined
% kFuncDeriv_mat: the cell array of the matlab functions of the derivatives of
% the kernel
% nb_samples: the number of samples taken randomly in the cube to estimate eta
% nb_angles: the number of angles taken deterministically to estimate eta
% in the convexity case
% bool_conv: if true, then etaDerivConv is computed
if nargin<6
    bool_conv=false;
end
centers=[0,0];

evalDeriv1=@(x,y) [kFuncDeriv_mat{2,1,2,1}(x,y,sigma)',...
    kFuncDeriv_mat{1,2,1,2}(x,y,sigma)'];

etaDeriv10=0; etaDeriv01=0;    etaDerivConv=0;
    loc_samples=(rand(nb_samples,2)-[1,1]/2).*deltaDesi+centers;
for j=1:nb_samples
    tmp_val=evalDeriv1(centers,centers)+evalDeriv1(loc_samples(j,:),loc_samples(j,:))...
    -evalDeriv1(centers,loc_samples(j,:))-evalDeriv1(loc_samples(j,:),centers);
    etaDeriv10=max(etaDeriv10,sqrt(tmp_val(1)));
    etaDeriv01=max(etaDeriv01,sqrt(tmp_val(2)));
end

if bool_conv
    angles=linspace(0,pi,nb_angles+1)'; angles=angles(1:end-1);
    u=[cos(angles),sin(angles)];

    evalDeriv2=@(x,y,r1,r2) [kFuncDeriv_mat{r1,r2,3,1}(x,y,sigma)',...
        kFuncDeriv_mat{r1,r2,2,2}(x,y,sigma)',...
        kFuncDeriv_mat{r1,r2,2,2}(x,y,sigma)',...
        kFuncDeriv_mat{r1,r2,1,3}(x,y,sigma)'];
    mat_DerivIdx={[3,1],[2,2];[2,2],[1,3]};

    loc_samples=(rand(nb_samples,2)-[1,1]/2).*deltaDesi+centers;
    for j=1:nb_samples
    for Ucount=1:nb_angles
        tmp_val=0;
        for Icount=1:2
        for Jcount=1:2
            valCell=mat_DerivIdx{Icount,Jcount};
            tmp_val=tmp_val+u(Ucount,Icount)*u(Ucount,Jcount)*...
                evalDeriv2(centers,centers,valCell(1),valCell(2))*...
                reshape(u(Ucount,:)'*u(Ucount,:),4,1)+u(Ucount,Icount)*u(Ucount,Jcount)*...
                evalDeriv2(loc_samples(j,:),loc_samples(j,:),valCell(1),valCell(2))*...
                reshape(u(Ucount,:)'*u(Ucount,:),4,1)-u(Ucount,Icount)*u(Ucount,Jcount)*...
                evalDeriv2(centers,loc_samples(j,:),valCell(1),valCell(2))*...
                reshape(u(Ucount,:)'*u(Ucount,:),4,1)-u(Ucount,Icount)*u(Ucount,Jcount)*...
                evalDeriv2(loc_samples(j,:),centers,valCell(1),valCell(2))*...
                reshape(u(Ucount,:)'*u(Ucount,:),4,1);
        end
        end
        etaDerivConv=max(etaDerivConv,sqrt(tmp_val));
    end
    end
end
end

%%
% clear all
% load('GaussianFuncDeriv.mat');
% deltaDesi=[0.219478155023285 0.233668941094316]/2;
% sigma=0.442;
% 
% nb_samples=100; nb_centers=1;
% centers=rand(nb_centers,2);
% 
% evalDeriv1=@(x,y) [kFuncDeriv_mat{2,1,2,1}(x,y,sigma)',...
%     kFuncDeriv_mat{1,2,1,2}(x,y,sigma)'];
% 
% etaDeriv10=zeros(nb_centers,1); etaDeriv01=zeros(nb_centers,1);
% for i=1:nb_centers
%     loc_samples=(rand(nb_samples,2)-[1,1]/2).*deltaDesi+centers(i,:);
% for j=1:nb_samples
%     tmp_val=evalDeriv1(centers(i,:),centers(i,:))+evalDeriv1(loc_samples(j,:),loc_samples(j,:))...
%     -evalDeriv1(centers(i,:),loc_samples(j,:))-evalDeriv1(loc_samples(j,:),centers(i,:));
%     etaDeriv10(i)=max(etaDeriv10(i),sqrt(tmp_val(1)));
%     etaDeriv01(i)=max(etaDeriv01(i),sqrt(tmp_val(2)));
% end
% end
% %
% nb_samples=50; nb_centers=1; nb_angles=30; 
% centers=rand(nb_centers,2); 
% angles=linspace(0,pi,nb_angles+1)'; angles=angles(1:end-1);
% u=[cos(angles),sin(angles)];
% 
% evalDeriv2=@(x,y,r1,r2) [kFuncDeriv_mat{r1,r2,3,1}(x,y,sigma)',...
%     kFuncDeriv_mat{r1,r2,2,2}(x,y,sigma)',...
%     kFuncDeriv_mat{r1,r2,2,2}(x,y,sigma)',...
%     kFuncDeriv_mat{r1,r2,1,3}(x,y,sigma)'];
% 
% etaDerivConv=zeros(nb_centers,1);
% mat_DerivIdx={[3,1],[2,2];[2,2],[1,3]};
% 
% for i=1:nb_centers
%     loc_samples=(rand(nb_samples,2)-[1,1]/2).*deltaDesi+centers(i,:);
% for j=1:nb_samples
% for Ucount=1:nb_angles
%     tmp_val=0;
%     for Icount=1:2
%     for Jcount=1:2
%         valCell=mat_DerivIdx{Icount,Jcount};
%         tmp_val=tmp_val+u(Ucount,Icount)*u(Ucount,Jcount)*...
%             evalDeriv2(centers(i,:),centers(i,:),valCell(1),valCell(2))*...
%             reshape(u(Ucount,:)'*u(Ucount,:),4,1)+u(Ucount,Icount)*u(Ucount,Jcount)*...
%             evalDeriv2(loc_samples(j,:),loc_samples(j,:),valCell(1),valCell(2))*...
%             reshape(u(Ucount,:)'*u(Ucount,:),4,1)-u(Ucount,Icount)*u(Ucount,Jcount)*...
%             evalDeriv2(centers(i,:),loc_samples(j,:),valCell(1),valCell(2))*...
%             reshape(u(Ucount,:)'*u(Ucount,:),4,1)-u(Ucount,Icount)*u(Ucount,Jcount)*...
%             evalDeriv2(loc_samples(j,:),centers(i,:),valCell(1),valCell(2))*...
%             reshape(u(Ucount,:)'*u(Ucount,:),4,1);
%         
%     end
%     end
%     etaDerivConv(i)=max(etaDerivConv(i),tmp_val);
% end
% end
% end