% The code below allow to reproduce the full numerical result for the 
% economic arm experiment of the article
% P.-C. Aubin-Frankowski, Z. Szabo "Handling Hard Affine SDP Shape
% Constraints in RKHSs ", JMLR 2022
% This partial experiment should take about 1min to run on a normal computer.

% CVX http://cvxr.com/cvx/ SHOULD BE INSTALLED, we used a free academic
% license for MOSEK, the solver used in this code 
% It can be retrieved at: https://www.mosek.com/downloads/

% CODE BELOW TO RUN ONCE SO THAT MOSEK+CVX WOULD WORK
% addpath C:\Users\pierr\Mosek\9.2\toolbox\R2015aom
% cd C:\Users\pierr\Desktop\cvx
% cvx_setup

%%Base set generation
clear all
load('dataLabourNorm-crop.mat');

tic
sigPrct=80;
sigma=sqrt(prctile(pdist(unique(data(:,1:2)),'squaredeuclidean'), sigPrct));%1.068675488891065;
X=data(:,1:2); Y=data(:,3); %Xadd=[0,0]; 

m1=min(data(:,1)); m2=min(data(:,2)); nbPts_per_dimCons=15;
deltaDesi=([2,2]-[m1,m2])/nbPts_per_dimCons;
xCoord=m1+deltaDesi(1)*(1/2+(0:nbPts_per_dimCons-1));
yCoord=m2+deltaDesi(2)*(1/2+(0:nbPts_per_dimCons-1));
[grid1add,grid2add] = meshgrid(xCoord,yCoord);
Xadd=[grid1add(:),grid2add(:)]; 

offsetGridmin=0.1;  offsetGridmax=2; nbPts_per_dimGrid=50;%max(X(:,1))+offsetGridmax
[grid1,grid2]=meshgrid(linspace(min(data(:,1))-offsetGridmin,4,nbPts_per_dimGrid),...
    linspace(min(data(:,2))-offsetGridmin,4,nbPts_per_dimGrid));
Xgrid=[grid1(:),grid2(:)]; 
ngrid=size(Xgrid,1); nadd=size(Xadd,1); nX=size(X,1); ntot=nX+nadd*6;

[etaDeriv10,etaDeriv01,etaDerivConv] = KRR_LabourEta(sigma,deltaDesi,kFuncDeriv_mat,50,20,true);
%% Computing large matrices from which to extract submatrices
[EvalXgridXtot, EvalD00XaddXtot, EvalD00XXtot, EvalD10XaddXtot, EvalD01XaddXtot,...
    EvalD20XaddXtot, EvalD11XaddXtot, EvalD11bisXaddXtot, EvalD02XaddXtot]...
    = KRR_LabourEvalFuncGaussian(sigma,kFuncDeriv_mat,X,Xadd,Xgrid);
tol=1E-7;
Gtot=[EvalD00XXtot; EvalD10XaddXtot; EvalD01XaddXtot;EvalD20XaddXtot;...
    EvalD11XaddXtot;EvalD11bisXaddXtot;EvalD02XaddXtot]+tol*eye(ntot);   
Gcons=chol(Gtot);
%% Cross-validation part (to repeat 20 times)

PropPts_CV=0.05; nb_folds=ceil(1/PropPts_CV);
nbPts_per_fold=floor(nX/nb_folds);
boolMonot=1; boolSOCMonot=0; boolConv=1; boolSOCConv=0;
factorM=2; factorC=20;%factorC=2 in first batch
for count=1:5
PropTrain=0.05; randIdx=randperm(nX); 
trainIdx=randIdx(1:floor(PropTrain*nX)); testIdx=randIdx(floor(PropTrain*nX)+1:end);
Xtrain=X(trainIdx,:); ntrain=length(trainIdx);
Xtest=X(testIdx,:); ntest=length(testIdx);
for boundN=[10]% 30 50 100 %[0.1 0.5 1 5 10 20 30 40 50 100 500 1E3] %

    tic
CVX_ploss=@(Xtrain,Xtest) CrossVal_KRR_LabourData_RMSE(Xtrain,Xtest,X,Xadd,Y,...
    etaDeriv10,etaDeriv01,etaDerivConv,...
    Gcons,EvalD00XXtot, EvalD10XaddXtot, EvalD01XaddXtot,...
    EvalD20XaddXtot, EvalD11XaddXtot, EvalD11bisXaddXtot, EvalD02XaddXtot,...
    factorM,factorC*boundN,boolMonot,boolSOCMonot,boolConv,boolSOCConv,boundN);          
results=CVX_ploss(Xtrain,Xtest);
% vals_cv = crossval(CVX_ploss,Xtrain,'k',nb_folds);%floor(PropPts_CV*nX)
%
saveFileName=['postCrossVal_table_results_Labour_wNorm_dummy.txt'];
fileID = fopen(saveFileName,'a');
fprintf(fileID,'%.4f \t %.4f \t %.4f \t %.2f \t %.2f \t %d \t %.2e \t %.2e \t %d  \t %d  \t %d  \t %d \t %d  \t %d \t %d \n',...%\t %.2e \t %d 
        results, PropPts_CV, sigPrct, boundN, boolMonot,boolSOCMonot,boolConv,boolSOCConv, factorM, factorC, count);
fclose(fileID);    
% for i=1:size(vals_cv,1)
% fprintf(fileID,'%.4f \t %.4f \t %.4f \t %.2f \t %.2f \t %d \t %.2e \t %.2e \t %d  \t %d  \t %d  \t %d \t %d  \t %d \n',...%\t %.2e \t %d 
%         vals_cv(i,:), PropPts_CV, sigPrct, boundN, boolMonot,boolSOCMonot,boolConv,boolSOCConv, factorM, factorC);
% end
% fclose(fileID);
elapsedTime=toc;
disp(['Finished the step for bound=' num2str(boundN) ' in ' num2str(elapsedTime) 's'])
end
end
disp(['Finished all the steps'])
% if savebool
%     saveFileName=['QR_Table_Results_2D.txt'];
%     fileID = fopen(saveFileName,'a');
%     fprintf(fileID,'%d \t %.2f \t %.2f \t %.2f \t %d \t %.2e \t %d \t %.2e \t %.2e \t %.2e\n',...
%             fileNb,test_ploss, best_sPrct, sigX, best_boundN, boundB);
%     fclose(fileID);
% end