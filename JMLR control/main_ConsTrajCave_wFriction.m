% The code below allow to reproduce the numerical result for the 
% control of submarine experiment of the article
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
close all
m=1; Asys=[0 1;0 -1/m]; Bsys=[0;1/m]; N=size(Asys,1); C=[1 0]'; nc=size(C,1);
%%
%  sig=0.1; nbRdmPts=20; nbDscPts=50; t_init=0; t_fin=1;
% scurr = rng;
% save('rng_Cave.mat','scurr');
load('rng_Cave.mat');
sig=0.08; nbRdmPts=20; nbDscPts=50; t_init=0; t_fin=1;
x1_init=0; x2_init=0; %x1_final=0; 
x_init=[x1_init x2_init]'; %x_inter=[.5 NaN NaN]'; xfinal=[0 NaN NaN]';
x_cons=x_init; listTcons=[t_init];
% nAdd= 20; listT=linspace(0,t_fin,nAdd); 

[listTbounds,xmin_arr,xmax_arr,gridTbis,valslow,valsup] = randomBorderGeneration(nbRdmPts,nbDscPts,sig,t_init, t_fin, scurr);%
% listTbounds=[0 0.1 0.2 0.3 0.4 0.6 1]; listXmin=[-1  -1 -1 1 2 0.8 0.8]; listXmax=[0.5 -0.7 1.5 2.2 3 2.5 2.5];
nbounds=length(listTbounds); 
% xmin_arr=-5*ones(nbounds,1); xmax_arr=5*ones(nbounds,1);
% for i=2:nbounds
%         xmin_arr(i)=max(listXmin(i-1),listXmin(i));
%         xmax_arr(i)=min(listXmax(i-1),listXmax(i));
% end

ntest=101; gridT=linspace(0,1,ntest); 

nAdd= length(listTbounds); listT=listTbounds;
listT=[listTcons,listT]; ntot=length(listT);
%% Gram matrices
[GX0,GX1]=GramianComputingVanLoanNonSym(listT,listT,Asys,Bsys);
[GXcons0,GXcons1]=GramianComputingVanLoanNonSym(listTbounds,listT,Asys,Bsys);
[GXtest0,GXtest1]=GramianComputingVanLoanNonSym(gridT,listT,Asys,Bsys);
%%
ratio=0;
GX=ratio*GX0+GX1;
tol=1E-10; GXsqrt=chol(GX+tol*eye(ntot*N));
GXcons=ratio*GXcons0+GXcons1;
GXtest=ratio*GXtest0+GXtest1;
%%
P_bool=[~isnan((x_cons));repmat(sum(logical(C),2),nAdd,1)];
P_bool_list_absent=find(~P_bool);
P_bool=repmat(sum(logical(C),2),nbounds,1);
P_bool_list_cons=find(P_bool);
%%
delta=(t_fin-t_init)/(2*nbDscPts);
eta_arr_cell = EtaComputingVanLoan(listTbounds, repmat(delta,1,nbounds), 2, Asys, Bsys, C, 0);
%%
cvx_begin
    cvx_solver mosek_2
    variables A(ntot*N,1)
    minimize(norm(GXsqrt*A))%
    subject to %
        0==A(P_bool_list_absent);
       xmin_arr+eta_arr_cell*norm(GXsqrt*A) <= GXcons(P_bool_list_cons,:)*A;
       xmax_arr-eta_arr_cell*norm(GXsqrt*A) >= GXcons(P_bool_list_cons,:)*A;
cvx_end
Abest=A;

eta=0;
cvx_begin
    cvx_solver mosek_2
    variables A(ntot*N,1)
    minimize(norm(GXsqrt*A))%
    subject to %
        0==A(P_bool_list_absent);
       xmin_arr+eta*norm(GXsqrt*A) <= GXcons(P_bool_list_cons,:)*A;
       xmax_arr-eta*norm(GXsqrt*A) >= GXcons(P_bool_list_cons,:)*A;
cvx_end

Abest0=A;%rand(ntot*N,1);

%% Moved constraint
Abest_alt=Abest;
Abest_alt(abs(Abest_alt)<1)=0;
Xtest=GXtest*Abest_alt; Xtest1=reshape(Xtest,2,[]);Xtest1=Xtest1(1,:)';
Xtest0eta=GXtest*Abest0; Xtest0eta1=reshape(Xtest0eta,2,[]);Xtest0eta1=Xtest0eta1(1,:)';

figure
hold on
plot(gridT,Xtest1,'r','LineWidth',2)
plot(gridT,Xtest0eta1,'--r','LineWidth',2)
plot([0],[x1_init],'ok','HandleVisibility','off')%,'Color',[0.6350 0.0780 0.1840]
quiver(0,0,0.1,0,'k','HandleVisibility','off')
% stairs(listTbounds'-delta,xmin_arr,'b','LineWidth',2)
% stairs(listTbounds'-delta,xmax_arr,'b','LineWidth',2,'HandleVisibility','off')
stairs(listTbounds'-delta,xmin_arr+eta_arr_cell*norm(GXsqrt*Abest),'g','LineWidth',2)
stairs(listTbounds'-delta,xmax_arr-eta_arr_cell*norm(GXsqrt*Abest),'g','LineWidth',2,'HandleVisibility','off')
plot(gridTbis,valslow,'k','LineWidth',1)
plot(gridTbis,valsup,'k','LineWidth',1,'HandleVisibility','off')
axis tight
%
%     
lgd=legend([{'Optimal trajectory with SOC constraints (ball covering)'};...
    {'Optimal trajectory with discretized constraints ($\eta=0$)'};...
%     {'Upper/Lower constraints $z_{low,m}$ and $z_{up,m}$'};...
    {'Upper/Lower constraints $z_{low,m}$ and $z_{up,m}\pm\eta_m \|\bar{f}\|_{K}$'};...
    {'Functions from Gaussian RKHS used to generate the bounds'}],'Interpreter','latex','FontSize',10);%Upper/Lower constraints $z_{low}$ and $z_{up}$
lgd.Location='southeast'; %
xlabel({'$t$'},'Interpreter','latex','FontSize',14); 
ylabel({'$z(t)$'},'Interpreter','latex','FontSize',14)
% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'Cave_TvsZ_Traj_movedConsJMLR','-dpdf')

%% Original constraint
Abest_alt=Abest;
Abest_alt(abs(Abest_alt)<1)=0;
Xtest=GXtest*Abest_alt; Xtest1=reshape(Xtest,2,[]);Xtest1=Xtest1(1,:)';
Xtest0eta=GXtest*Abest0; Xtest0eta1=reshape(Xtest0eta,2,[]);Xtest0eta1=Xtest0eta1(1,:)';

figure
hold on
plot(gridT,Xtest1,'r','LineWidth',2)
plot(gridT,Xtest0eta1,'--r','LineWidth',2)
plot([0],[x1_init],'ok','HandleVisibility','off')%,'Color',[0.6350 0.0780 0.1840]
quiver(0,0,0.1,0,'k','HandleVisibility','off')
stairs(listTbounds'-delta,xmin_arr,'b','LineWidth',2)
stairs(listTbounds'-delta,xmax_arr,'b','LineWidth',2,'HandleVisibility','off')
% stairs(listTbounds'-delta,xmin_arr+eta_arr_cell*norm(GXsqrt*Abest),'g','LineWidth',2)
% stairs(listTbounds'-delta,xmax_arr-eta_arr_cell*norm(GXsqrt*Abest),'g','LineWidth',2,'HandleVisibility','off')
% plot(gridTbis,valslow,'k','LineWidth',1)
% plot(gridTbis,valsup,'k','LineWidth',1,'HandleVisibility','off')
plot([0],[x1_init],'ok','HandleVisibility','off')%,'Color',[0.6350 0.0780 0.1840]
quiver(0,0,0.1,0,'k','HandleVisibility','off')
axis tight

lgd=legend([...
    {'Optimal trajectory with SOC constraints (ball covering)'};...
    {'Optimal trajectory with discretized constraints ($\eta=0$)'};...
    {'Upper/Lower constraints $z_{low,m}$ and $z_{up,m}$'};...
%     {'Upper/Lower constraints $z_{low,m}$ and $z_{up,m}\pm\eta_m \|\bar{x}(\cdot)\|_{K}$'};...
%     {'Functions from Gaussian RKHS used to generate the bounds'}
    ],'Interpreter','latex','FontSize',10);%Upper/Lower constraints $z_{low}$ and $z_{up}$
lgd.Location='southeast'; %
xlabel({'$t$'},'Interpreter','latex','FontSize',14); 
ylabel({'$z(t)$'},'Interpreter','latex','FontSize',14)
% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'Cave_TvsZ','-dpdf')
% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'Cave_TvsZ_Traj_ConsJMLR','-dpdf')