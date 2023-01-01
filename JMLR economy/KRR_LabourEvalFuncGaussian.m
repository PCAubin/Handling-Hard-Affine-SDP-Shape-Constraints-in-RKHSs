function [EvalXgridXtot, EvalD00XaddXtot, EvalD00XXtot, EvalD10XaddXtot, EvalD01XaddXtot,...
    EvalD20XaddXtot, EvalD11XaddXtot, EvalD11bisXaddXtot, EvalD02XaddXtot]...
    = KRR_LabourEvalFuncGaussian(sigma,kFuncDeriv_mat,X,Xadd,Xgrid)
%KRR_LABOUREVALFUNCGAUSSIAN computes the various Gram matrices of the kernel and its derivatives appearing in
%the problem based on a kernel of bandwith sigma
% EvalD20XaddXtot for instance corresponds to the second derivative w.r.t.
% the first argument of the kernel between Xadd and Xtot (the total number
% of points, i.e. the measurement points X and the constraints Xadd counted
% six times)
ngrid=size(Xgrid,1); nadd=size(Xadd,1); n=size(X,1); ntot=n+nadd*6;

evalDeriv=@(x,r1,r2) [kFuncDeriv_mat{r1,r2,1,1}(x,X,sigma)',...
    kFuncDeriv_mat{r1,r2,2,1}(x,Xadd,sigma)',...
    kFuncDeriv_mat{r1,r2,1,2}(x,Xadd,sigma)',...
    kFuncDeriv_mat{r1,r2,3,1}(x,Xadd,sigma)',...
    kFuncDeriv_mat{r1,r2,2,2}(x,Xadd,sigma)',...
    kFuncDeriv_mat{r1,r2,2,2}(x,Xadd,sigma)',...
    kFuncDeriv_mat{r1,r2,1,3}(x,Xadd,sigma)'];

EvalD00XXtot=zeros(n,ntot);
EvalD00XaddXtot=zeros(nadd,ntot);
EvalD10XaddXtot=zeros(nadd,ntot); 
EvalD01XaddXtot=zeros(nadd,ntot);
EvalD20XaddXtot=zeros(nadd,ntot);
EvalD11XaddXtot=zeros(nadd,ntot);
EvalD11bisXaddXtot=zeros(nadd,ntot);
EvalD02XaddXtot=zeros(nadd,ntot);
EvalXgridXtot=zeros(ngrid,ntot);
for i=1:nadd
    EvalD00XaddXtot(i,:)=evalDeriv(Xadd(i,:),1,1);
    EvalD10XaddXtot(i,:)=evalDeriv(Xadd(i,:),2,1);
    EvalD01XaddXtot(i,:)=evalDeriv(Xadd(i,:),1,2);
    EvalD20XaddXtot(i,:)=evalDeriv(Xadd(i,:),3,1);
    EvalD11XaddXtot(i,:)=evalDeriv(Xadd(i,:),2,2);
    EvalD11bisXaddXtot(i,:)=evalDeriv(Xadd(i,:),2,2);
    EvalD02XaddXtot(i,:)=evalDeriv(Xadd(i,:),1,3);
end
for i=1:n
    EvalD00XXtot(i,:)=evalDeriv(X(i,:),1,1);
end
for i=1:ngrid
     EvalXgridXtot(i,:)=evalDeriv(Xgrid(i,:),1,1);
end
end

