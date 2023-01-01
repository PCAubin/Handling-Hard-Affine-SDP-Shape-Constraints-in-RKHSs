function [newTDeltaEta_mat] = ComputeCovering_Delta(TDeltaEta_mat,gamma,lbda,hlap)
%Only for translation-invariant kernels, the function outputs the eta
%quantities newTDeltaEta_mat, for a given Matern kernel with bandwidth lbda, 
%after one step of grid refinement based on the previous values
%TDeltaEta_mat. gamma is the shrinking factor. This function uses the fact
%that in 1D and for a Matérn kernel all the computations are explicit.

newEta_arr=gamma*TDeltaEta_mat(:,3);
newDeltaMax_arr=-log(1-newEta_arr.^2/2)/lbda;%%formula for Matérn kernel
nbptsToChange=size(TDeltaEta_mat,1);
newDelta_arr=TDeltaEta_mat(:,2)./ceil(TDeltaEta_mat(:,2)./newDeltaMax_arr);
newTDeltaEta_cellmat=cell(nbptsToChange,1);
for i=1:nbptsToChange
tempTDE=TDeltaEta_mat(i,1)-TDeltaEta_mat(i,2)+newDelta_arr(i):2*newDelta_arr(i):TDeltaEta_mat(i,1)+TDeltaEta_mat(i,2)-newDelta_arr(i);
newTDeltaEta_cellmat{i}=[tempTDE',repmat(newDelta_arr(i),length(tempTDE),1),repmat(sqrt(2*(1-hlap(newDelta_arr(i),lbda))),length(tempTDE),1)];
end
newTDeltaEta_mat=cell2mat(newTDeltaEta_cellmat);
end

% newEta_arr=gamma*TDeltaEta_mat(:,3);
% newDelta_arr=-log(1-newEta_arr.^2/2)/lbda;%%formula for Matérn kernel
% nbpts=size(TDeltaEta_mat,1);
% newTDeltaEta_cellmat=cell(nbpts,1);
% for i=1:nbpts
% tempTDE=TDeltaEta_mat(i,1)-TDeltaEta_mat(i,2)+newDelta_arr(i):newDelta_arr(i):TDeltaEta_mat(i,1)+TDeltaEta_mat(i,2);
% newTDeltaEta_cellmat{i}=[tempTDE',repmat(newDelta_arr(i),length(tempTDE),1),repmat(newEta_arr(i),length(tempTDE),1)];
% end
% newTDeltaEta_mat=cell2mat(newTDeltaEta_cellmat);