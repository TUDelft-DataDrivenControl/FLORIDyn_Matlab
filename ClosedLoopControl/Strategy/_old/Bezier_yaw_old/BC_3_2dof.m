function yaw_off_setpoints = BC_3_2dof(yaw_rate_limit_norm,yaw_n)
%BC_3_2DOF Cubic Bézier-Curve 
% INPUTS
% yaw_rate_limit_norm = [double] yaw rate limit multiplied by the time
%                               window. Has to be larger than 0.
% yaw_n = [n x 2] normalized choices for the yaw angles \in [0,1] 
% OUTPUTS
% yaw_off_setpoints = [n x 2] yaw offset setpoints for the near future.
%% Convert [0,1] to yaw setpoints
yaw_off_setpoints = zeros(size(yaw_n));

yaw_off_setpoints(:,1) = (yaw_n(:,1) * 2 - 1) * yaw_rate_limit_norm/3;

h1 = -1/3*( -sqrt( yaw_rate_limit_norm^2 - ...
    3*yaw_rate_limit_norm * yaw_off_setpoints(:,1))-...
    yaw_rate_limit_norm - 3 * yaw_off_setpoints(:,1));
h2 = 1/3*( -sqrt( yaw_rate_limit_norm^2 + ...
    3*yaw_rate_limit_norm * yaw_off_setpoints(:,1))-...
    yaw_rate_limit_norm + 3 * yaw_off_setpoints(:,1));

yaw_off_setpoints(:,2) = h2 + yaw_n(:,2).*(h1 - h2);
end

