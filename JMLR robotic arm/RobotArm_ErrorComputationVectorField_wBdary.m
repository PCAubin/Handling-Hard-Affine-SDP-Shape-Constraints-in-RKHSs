function error_arr = RobotArm_ErrorComputationVectorField_wBdary(v_true,v_approx,vals_bdry)
%ERRORCOMPUTATIONVECTORFIELD outputs an array (4*nb_points) corresponding
%to four classical metrics to assess the error between a true vector field
%v_true (since we have vector outputs) and its approximation v_approx. The
%L_2 and L_infty error are followed by two metrics to assess violation of
%constraints which are here requiring nonnegativity of some quantity
%vals_bdry. For the violation, we compute L_1 and L_infty errors. The Q2
%error is similar the L2 error over the fields but with some extra normalization. In the
%article only L2Err and BdryErrL1 are reported.

nb_points=size(v_true,1);

L2Err=sum(vecnorm(v_true-v_approx,2,2).^2)/nb_points;
LinftyErr=max(vecnorm(v_true-v_approx,2,2));
BdryErrLinfty=-min(min(vals_bdry),0);
BdryErrL1=-sum(min(vals_bdry,0))/size(vals_bdry,1);

Q2=1-sum(vecnorm(v_true-v_approx,2,2).^2)/sum(vecnorm(v_true-mean(v_true),2,2).^2);
error_arr=[Q2,L2Err,LinftyErr,BdryErrLinfty,BdryErrL1];
end

