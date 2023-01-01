function output=robot_arm(x)
%ROBOT_ARM gives the "tool tip position", 2D-position and 1D-orientation of
%the tip of a robot arm in the state described by x (which combines lengths
%of links and angles of joints).
    m = size(x,2)/2;
%     assert n % 2 == 0, 'only even number input for robot arm function'
    L = x(:,1:m);
    theta = 2*pi*x(:,m+1:end);

    theta_sums = cumsum(theta,2);
    y=[sum(L.*cos(theta_sums),2),sum(L.*sin(theta_sums),2)];
    z=sin(sum(theta,2));
    output=[y,z];
end
%https://github.com/cagrell/gp_constr/blob/master/Example_2.ipynb