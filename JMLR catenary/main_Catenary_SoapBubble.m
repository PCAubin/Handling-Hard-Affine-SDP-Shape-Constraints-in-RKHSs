% The code below allow to reproduce the numerical result for the 
% catenary experiment of the article
% P.-C. Aubin-Frankowski, Z. Szabo "Handling Hard Affine SDP Shape
% Constraints in RKHSs ", JMLR 2022
% This experiment should take about 1 min to run on a normal computer.
% The code will produce a table (to be saved in .mat) storing all the
% relevant quantities which can then be plotted in plot_Catenary_SoapBubble.

% CVX http://cvxr.com/cvx/ SHOULD BE INSTALLED, we used a free academic
% license for MOSEK, the solver used in this code 
% It can be retrieved at: https://www.mosek.com/downloads/

% CODE BELOW TO RUN ONCE SO THAT MOSEK+CVX WOULD WORK
% addpath C:\Users\pierr\Mosek\9.2\toolbox\R2015aom
% cd C:\Users\pierr\Desktop\cvx
% cvx_setup

clear all 

ntest=400+3; gridT=linspace(0,1,ntest);
hlap = @(u,lbda) exp(-lbda * abs(u)); h=hlap; lbda=5;
xinit=0; xend=xinit; xmid=1.5; xcons=0.5; tconsInit=0.2; tconsEnd=0.8;
DltaInit=0.01; ListTconsInit=unique([tconsInit+DltaInit:2*DltaInit:tconsEnd-DltaInit,tconsEnd-DltaInit]);
ncons=length(ListTconsInit); etaInit=sqrt(2*(1-hlap(DltaInit,lbda)));
TDeltaEta_matInit=[ListTconsInit',repmat(DltaInit,ncons,1),repmat(etaInit,ncons,1)]; 
listTcons=[0 0.5 1];

listTInit=[listTcons,ListTconsInit]; nInit=length(listTInit);

%%
nbIter=5;gamma=0.8;
% nbIter=30; % parameter chosen in the article, takes 1h to run.
TDeltaEta_cellmat=cell(nbIter+1,6); % TED, listT, A, distance to cons, consIdx, cputime, optval
TDeltaEta_cellmat{1,1}=TDeltaEta_matInit;
TDeltaEta_cellmat{1,2}=listTInit;

tolCons=1E-8;

for j=2:nbIter+1 
    n=length(TDeltaEta_cellmat{j-1,2}); ncons=length(TDeltaEta_cellmat{j-1,1}(:,1));
    if (n>500)||(n==0)
        break;
    else
    GX=h(repmat(TDeltaEta_cellmat{j-1,2},n,1)-repmat(TDeltaEta_cellmat{j-1,2}',1,n),lbda); 
    GXsqrt=chol(GX);  
cvx_begin quiet
    cvx_solver mosek
    variables A(n,1)
    minimize(norm(GXsqrt*A))
    subject to 
        xinit==GX(1,:)*A;
        xmid==GX(2,:)*A;
        xend==GX(3,:)*A; 
        for i=length(listTcons)+1:n
        TDeltaEta_cellmat{j-1,1}(i-length(listTcons),3)*norm(GXsqrt*A)+xcons <= GX(i,:)*A;
        end
cvx_end
% Abest=A; Xtest=GXtest*Abest; 
TDeltaEta_cellmat{j-1,3}=A;
TDeltaEta_cellmat{j-1,4}=TDeltaEta_cellmat{j-1,1}(:,3)*norm(GXsqrt*A)+xcons - GX(length(listTcons)+1:n,:)*A;
TDeltaEta_cellmat{j-1,6}=cvx_cputime;
TDeltaEta_cellmat{j-1,7}=cvx_optval;

consIdx=find(abs(TDeltaEta_cellmat{j-1,1}(:,3)*norm(GXsqrt*A)+xcons - GX(length(listTcons)+1:n,:)*A)<tolCons);
TDeltaEta_cellmat{j-1,5}=TDeltaEta_cellmat{j-1,1}(consIdx,1:2);

TDeltaEta_cellmat{j,1}=[TDeltaEta_cellmat{j-1,1};...
    ComputeCovering_Delta(TDeltaEta_cellmat{j-1,1}(consIdx,:),gamma,lbda,hlap)];
TDeltaEta_cellmat{j,1}(consIdx,:)=[];
TDeltaEta_cellmat{j,2}=[listTcons,TDeltaEta_cellmat{j,1}(:,1)'];
    end
end
disp("End computation ball covering")
save(['Catenary_G' num2str(gamma*10) '_L' num2str(lbda) '_I' num2str(nbIter-1) '_balls.mat'], 'TDeltaEta_cellmat')

%% Hyperplanes
nbIter=5;gamma=0.8;
% nbIter=40; % parameter chosen in the article, takes 1h to run.
TDeltaEta_cellmat=cell(nbIter+1,6); % TED, listT, A, distance to cons, conxIdx cputime, optval
TDeltaEta_cellmat{1,1}=TDeltaEta_matInit;
TDeltaEta_cellmat{1,2}=listTInit;
TDeltaEta_cellmat{1,4}=hlap(TDeltaEta_matInit(:,2),lbda);
tolCons=1E-7;

for j=2:nbIter+1 
    n=length(TDeltaEta_cellmat{j-1,2}); ncons=length(TDeltaEta_cellmat{j-1,1}(:,1));
    if n>500
        break;
    else
    GX=h(repmat(TDeltaEta_cellmat{j-1,2},n,1)-repmat(TDeltaEta_cellmat{j-1,2}',1,n),lbda); 
    GXsqrt=chol(GX);  
cvx_begin quiet
    cvx_solver mosek
    variables A(n,1) B(n,1)
    minimize(norm(GXsqrt*A))
    subject to 
        xinit==GX(1,:)*A;
        xmid==GX(2,:)*A;
        xend==GX(3,:)*A; 
        B>=0;
        for i=length(listTcons)+1:n
        boolVec=zeros(n,1); boolVec(i)=1;
        norm(GXsqrt*(A-B(i)*boolVec))+xcons <= B(i)*TDeltaEta_cellmat{j-1,4}(i-length(listTcons));
        end
cvx_end
TDeltaEta_cellmat{j-1,3}=A;
TDeltaEta_cellmat{j-1,6}=cvx_cputime;
TDeltaEta_cellmat{j-1,7}=cvx_optval;

listValsConsHyper=zeros(n-length(listTcons),1);
for i=length(listTcons)+1:n
    boolVec=zeros(n,1); boolVec(i)=1;
    listValsConsHyper(i-length(listTcons))=norm(GXsqrt*(A-B(i)*boolVec))+xcons - B(i)*TDeltaEta_cellmat{j-1,4}(i-length(listTcons));
end
        
consIdx=find(abs(listValsConsHyper)<tolCons);
TDeltaEta_cellmat{j-1,5}=TDeltaEta_cellmat{j-1,1}(consIdx,1:2);

TDeltaEta_cellmat{j,1}=[TDeltaEta_cellmat{j-1,1};...
    ComputeCovering_Delta(TDeltaEta_cellmat{j-1,1}(consIdx,:),gamma,lbda,hlap)];
TDeltaEta_cellmat{j,1}(consIdx,:)=[];
TDeltaEta_cellmat{j,2}=[listTcons,TDeltaEta_cellmat{j,1}(:,1)'];
TDeltaEta_cellmat{j,4}=hlap(TDeltaEta_cellmat{j,1}(:,2),lbda);
    end
end
%optval=1.5446
disp("End computation hyperplanes")
TDeltaEta_cellmat_hyper=TDeltaEta_cellmat;
% save(['Catenary_G' num2str(gamma*10) '_L' num2str(lbda) '_I' num2str(nbIter-1) '_halfspaces.mat'], 'TDeltaEta_cellmat')