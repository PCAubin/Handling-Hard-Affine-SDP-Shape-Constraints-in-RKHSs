function Sig_corelMat_th=RobotArm_estimateCovOutput(dim_in)
% RobotArm_estimateCovOutput outputs a matrix Sig_corelMat_th (3x3) which 
% is the empirical covariance of the 3D outputs of robot_arm, estimated over
% 1000 points taken at random through lhsdesign.
x_design=lhsdesign(1000,dim_in,'Criterion','maximin','Iterations',1000);
y_design =robot_arm(x_design);
Sig_corelMat_th=cov(y_design);
end