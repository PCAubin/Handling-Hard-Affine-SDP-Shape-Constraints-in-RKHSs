% The code below allow to reproduce a minimal numerical result for the 
% economic arm experiment of the article
% P.-C. Aubin-Frankowski, Z. Szabo "Handling Hard Affine SDP Shape
% Constraints in RKHSs ", JMLR 2022
% This partial experiment should take about 20s to run on a normal computer.
% The code will produce one of the 3D plots of the article depending on the
% constraints one chose to be active.

% CVX http://cvxr.com/cvx/ SHOULD BE INSTALLED, we used a free academic
% license for MOSEK, the solver used in this code 
% It can be retrieved at: https://www.mosek.com/downloads/

% CODE BELOW TO RUN ONCE SO THAT MOSEK+CVX WOULD WORK
% addpath C:\Users\pierr\Mosek\9.2\toolbox\R2015aom
% cd C:\Users\pierr\Desktop\cvx
% cvx_setup

clear all
load('dataLabourNorm-crop.mat');
load('GaussianFuncDeriv.mat');
tic
sigma=0.442;%;
% sigma=1.0705;%75
sigma=1.068675488891065;%sqrt(prctile(pdist(unique(data(:,1:2)),'squaredeuclidean'), 80));%1.3147;%80th decile of full data
% sigma=2.2975;%90
nbPtskept=20; idxPts=1:nbPtskept;%randperm(size(data,1),nbPtskept);

X=data(idxPts,1:2); Y=data(idxPts,3); %Xadd=[0,0]; 

m1=min(data(:,1)); m2=min(data(:,2)); nbPts_per_dimCons=15;

deltaDesi=([2,2]-min(X))/nbPts_per_dimCons;
xCoord=m1+deltaDesi(1)*(1/2+(0:nbPts_per_dimCons-1));
yCoord=m2+deltaDesi(2)*(1/2+(0:nbPts_per_dimCons-1));

[grid1add,grid2add] = meshgrid(xCoord,yCoord);
% [grid1add,grid2add]=meshgrid(linspace(min(X(:,1)),max(X(:,1)),10),linspace(min(X(:,2)),max(X(:,2)),10));
Xadd=[grid1add(:),grid2add(:)]; 

offsetGridmin=0.1;  offsetGridmax=2; nbPts_per_dimGrid=50;%max(X(:,1))+offsetGridmax
[grid1,grid2]=meshgrid(linspace(min(X(:,1))-offsetGridmin,4,nbPts_per_dimGrid),...
    linspace(min(X(:,2))-offsetGridmin,4,nbPts_per_dimGrid));
Xgrid=[grid1(:),grid2(:)]; 
ngrid=size(Xgrid,1); nadd=size(Xadd,1); n=size(X,1); ntot=n+nadd*6;

[etaDeriv10,etaDeriv01,etaDerivConv] = KRR_LabourEta(sigma,deltaDesi,kFuncDeriv_mat,50,20,true);
% 10 pts [1.1,6.3] 20 pts [0.56,3.4] 30 pts [0.36,2.3] 40 pts [0.28,1.7]
%%
[EvalXgridXtot, EvalD00XaddXtot, EvalD00XXtot, EvalD10XaddXtot, EvalD01XaddXtot,...
    EvalD20XaddXtot, EvalD11XaddXtot, EvalD11bisXaddXtot, EvalD02XaddXtot]...
    = KRR_LabourEvalFuncGaussian(sigma,kFuncDeriv_mat,X,Xadd,Xgrid);
%%
tol=1E-7;
Gtot=[EvalD00XXtot; EvalD10XaddXtot; EvalD01XaddXtot;EvalD20XaddXtot;...
    EvalD11XaddXtot;EvalD11bisXaddXtot;EvalD02XaddXtot]+tol*eye(ntot);   
Gcons=chol(Gtot);
%%
tic
normBound=2E1; 
etaMonot=1E-2; %1
etaConv=3E-2; %40
factorM=2; factorC=4E1;
%CHANGE BELOW TO (DE)ACTIVATE A CONSTRAINT
boolMonot=1; boolSOCMonot=1; boolConv=1; boolSOCConv=1;

if boolConv
    cvx_begin sdp
elseif ~boolConv
    cvx_begin
end
    cvx_precision low
    cvx_solver mosek_2
    variables A(ntot,1)
    minimize(norm(Y-EvalD00XXtot*A))
    subject to
       norm(Gcons*A)<=normBound;
       if boolMonot
       EvalD10XaddXtot*A <= -boolSOCMonot*etaDeriv10*norm(Gcons*A)/factorM;
       EvalD01XaddXtot*A <= -boolSOCMonot*etaDeriv01*norm(Gcons*A)/factorM;
       elseif ~boolMonot
           A(n+1:n+2*nadd)==0;
       end           
       if boolConv
       for i=1:nadd
%            [EvalD20XaddXtot(i,:)*A EvalD11XaddXtot(i,:)*A;...
%                EvalD11bisXaddXtot(i,:)*A EvalD02XaddXtot(i,:)*A] - boolSOCConv*etaConv*eye(2) == semidefinite(2);
           [EvalD20XaddXtot(i,:)*A EvalD11XaddXtot(i,:)*A;...
               EvalD11bisXaddXtot(i,:)*A EvalD02XaddXtot(i,:)*A]>=boolSOCConv*etaDerivConv*normBound*eye(2)/factorC;
       end
       elseif ~boolConv
           A(n+2*nadd+1:end)==0;
       end
cvx_end
%
vals=EvalXgridXtot*A;
figure
hold on
grid on
ZL = [-3,3];
caxis(ZL)
scatter3(X(:,1),X(:,2),ZL(1)*ones(n,1),'ok')
scatter3(Xadd(:,1),Xadd(:,2),ZL(1)*ones(nadd,1),'or')
scatter3(X(:,1),X(:,2),Y,'ok','filled')
scatter3(data(:,1),data(:,2),data(:,3),repmat(5,size(data,1),1),'ok','filled')
scatter3(Xadd(:,1),Xadd(:,2),EvalD00XaddXtot*A,'or','filled')
for i=1:n
plot3([X(i,1),X(i,1)],[X(i,2),X(i,2)],[ZL(1) Y(i)],'k','LineWidth',1)
% plot3(X(i,1),X(i,2),ZL(1),'bo')
end

surf(grid1,grid2,reshape(vals,size(grid1,1),size(grid1,2)),'FaceAlpha',.7)
% alpha(surfPlt,0.7);

xlim([min(Xgrid(:,1)) max(Xgrid(:,1))])
ylim([min(Xgrid(:,2)) max(Xgrid(:,2))])
zlim(ZL)

xlabel('Capital')
ylabel('Labour')
zlabel('-log(Output)')
view([50 30])

% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,['KRR_LabourData_' num2str(n) 'Pts_wMCCons'],'-dpdf')

% saveas(gcf,['LabourData_Monot' num2str(boolMonot) '_SOC' num2str(boolSOCMonot)...
%      '_Conv' num2str(boolConv) '_SOC' num2str(boolSOCConv) '.png'])
elapsedTime=toc;