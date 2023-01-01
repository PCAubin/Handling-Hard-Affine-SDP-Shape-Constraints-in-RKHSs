function [BigMat0,BigMat1] = GramianComputingVanLoanNonSym(listT1,listT2,A,B)
% GRAMIANCOMPUTINGVANLOAN outputs two large matrices that are the concatenation of
% the Gramians K(listT1_i,listT2_j) of lisT1*listT2, the output is obtained
% through a "cell2mat" of nn1*nn2 cell_array of N*N matrices.
% BigMat0: (nn1xN)*(nn2xN) for the uncontrolled part K_0 of the kernel
% BigMat1: (nn1xN)*(nn2xN) for the controlled part K_1 of the kernel
% The matrices K(s,t_m) are approximated using Van Loan's trick described in "Computing Integrals Involving the
% Matrix Exponential, CHARLES F. VAN LOAN, TAC,1978". We use the same notations F2 and G2 as in the paper.
% Inputs:
% listT1: 1xnn1, list of time points where the kernel matrices are evaluated
% listT2: 1xnn2, list of time points where the kernel matrices are centered
% A: NxN, B:NxP matrices defining x'=Ax+Bu
% NOTE THAT THE GRAM MATRICES ARE COMPUTED for R=eye(P), Q=0 with the LQR
% classical notations x'Qx+u'Ru.
% Assumes that the time lists are ordered.


tic
N=size(A,1); nn1=length(listT1); nn2=length(listT2);
BigMat_cellArr0=cell(nn1,nn1); BigMat_cellArr1=cell(nn1,nn1);
BigMat_cellArrF2_T1=cell(nn1,1); BigMat_cellArrG2_T1=cell(nn1,1);
BigMat_cellArrF2_T2=cell(nn2,1); BigMat_cellArrG2_T2=cell(nn2,1);

for i=1:nn1 %Van Loan's technique for computing Gramians
        temp = expm([A B*B';zeros(N,N) -A']*listT1(i));
        BigMat_cellArrF2_T1{i}=temp(1:N,1:N);
        BigMat_cellArrG2_T1{i}=temp(1:N,N+1:2*N);
end
for i=1:nn2 %Van Loan's technique for computing Gramians
        temp = expm([A B*B';zeros(N,N) -A']*listT2(i));
        BigMat_cellArrF2_T2{i}=temp(1:N,1:N);
        BigMat_cellArrG2_T2{i}=temp(1:N,N+1:2*N);
end
for i=1:nn1 %Computing upperhalf with Van Loan's technique for computing Gramians
    for j=1:nn2
        BigMat_cellArr0{i,j}=BigMat_cellArrF2_T1{i}*BigMat_cellArrF2_T2{j}';
        if listT1(i)<= listT2(j)
            BigMat_cellArr1{i,j}=BigMat_cellArrG2_T1{i}*BigMat_cellArrF2_T2{j}';
        else
            BigMat_cellArr1{i,j}=BigMat_cellArrF2_T1{i}*BigMat_cellArrG2_T2{j}';
        end
    end
end

elapsedTime=toc;
disp(['Van Loan: finished Gram ' num2str(elapsedTime) 's']);
BigMat0 = cell2mat(BigMat_cellArr0);
BigMat1 = cell2mat(BigMat_cellArr1);

end

