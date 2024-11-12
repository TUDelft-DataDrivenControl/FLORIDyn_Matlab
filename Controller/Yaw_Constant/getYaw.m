function yaw = getYaw(ConYawData,iT,~)
%GETYAW Return the yaw angle for the requested turbine(s)
% ======================================================================= %
% Individual Turbine value implementation
%   requires a .csv in the simulation folder called Control_YawInterp.csv
%   where each row is a
%       time, yaw_T0, yaw_T1, ... yaw_Tn
%   setpoint in time. The values are interploated linearly between the
%   setpoints.
% ======= Input ======
% ConYawData = (t,yaw_T0, yaw_T1, ... yaw_Tn)
% iT         = Index/Indeces of the turbines
% t          = time of request
% ======================================================================= %

yaw = ones(size(iT'))*ConYawData(1,1);
end

