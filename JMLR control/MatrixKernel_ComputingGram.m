function BigMat = MatrixKernel_ComputingGram(X1,X2,Kfunc)
% MatrixKernel_ComputingGram 
% X1 of size nn1*dim, X2 of size nn1*dim
tic
nn1=size(X1,1); nn2=length(X2,1); dim=size(X1,2);
BigMat=zeros(nn1*dim,nn2*dim);
for i=1:nn1
    for j=1:nn2
        BigMat(1+(i-1)*dim:i*dim,1+(j-1)*dim:j*dim)=Kfunc(X1(i,:)',X2(j,:)');
    end
end
elapsedTime=toc;
disp(['Finished Gram ' num2str(elapsedTime) 's']);
end

