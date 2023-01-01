function result = CrossVal_KRR_LabourData_RMSE(Xtrain,Xtest,X,Xadd,Y,...
    etaDeriv10,etaDeriv01,etaDerivConv,...
    Gcons,EvalD00XXtot, EvalD10XaddXtot, EvalD01XaddXtot,...
    EvalD20XaddXtot, EvalD11XaddXtot, EvalD11bisXaddXtot, EvalD02XaddXtot,...
    factorM,factorC,boolMonot,boolSOCMonot,boolConv,boolSOCConv,boundN)
%CROSSVAL_KRR_LABOURDATA_RMSE outputs the results table to store in a .txt
%for each of the cross-validation experiments. It extracts from the Gram
%matrices EvalD00XXtot etc the proper submatrices to then use CVX to
%produce the solution. The various boolMonot etc correspond to the
%constraints that are active or not.

trainIdx = ismember(X,Xtrain,'rows');
testIdx = ismember(X,Xtest,'rows');
testIdx_notTrain=testIdx-trainIdx; testIdx_notTrain(testIdx_notTrain<0)=0;
testIdx_notTrain=logical(testIdx_notTrain);

ntrain=sum(trainIdx); ntest=sum(testIdx); 
nadd=size(Xadd,1); nX=size(X,1); ntot=nX+nadd*6;

compTestIdx=ones(nX,1)-trainIdx-testIdx; compTestIdx(compTestIdx<0)=0;
compTestIdx=logical(compTestIdx);
ncomptest=sum(compTestIdx);
%ntest+ntrain+ncomptest-nX
% Gtot=[EvalD00XXtot; EvalD10XaddXtot; EvalD01XaddXtot;EvalD20XaddXtot;...
%     EvalD11XaddXtot;EvalD11bisXaddXtot;EvalD02XaddXtot]+tol*eye(ntot);   

if boolConv
    cvx_begin sdp quiet
elseif ~boolConv
    cvx_begin quiet
end
    cvx_precision low
    cvx_solver mosek_2
    variables A(ntot,1)
    minimize(norm(Y(trainIdx,:)-EvalD00XXtot(trainIdx,:)*A))
    subject to
        A(testIdx_notTrain)==0;
       norm(Gcons*A)<=boundN;
       if boolMonot
       EvalD10XaddXtot*A <= -boolSOCMonot*etaDeriv10*norm(Gcons*A)/factorM;
       EvalD01XaddXtot*A <= -boolSOCMonot*etaDeriv01*norm(Gcons*A)/factorM;
       elseif ~boolMonot
           A(nX+1:nX+2*nadd)==0;
       end           
       if boolConv
       for i=1:nadd
           [EvalD20XaddXtot(i,:)*A EvalD11XaddXtot(i,:)*A;...
               EvalD11bisXaddXtot(i,:)*A EvalD02XaddXtot(i,:)*A]>=boolSOCConv*etaDerivConv*boundN*eye(2)/factorC;
       end
       elseif ~boolConv
           A(nX+2*nadd+1:end)==0;
       end
cvx_end

errorYTrain=norm(Y(trainIdx,:)-EvalD00XXtot(trainIdx,:)*A)^2/ntrain;
errorYTest=norm(Y(testIdx,:)-EvalD00XXtot(testIdx,:)*A)^2/ntest;
errorYcompTest=norm(Y(compTestIdx,:)-EvalD00XXtot(compTestIdx,:)*A)^2/ncomptest;
result=[errorYTest, errorYcompTest, errorYTrain, cvx_cputime, norm(Gcons*A)];%, cvx_slvitr, cvx_slvtol

% display("One cross-validation step done")
end