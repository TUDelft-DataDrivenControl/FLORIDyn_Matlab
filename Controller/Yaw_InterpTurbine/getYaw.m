function yaw = getYaw(ConYawData,iT,t)
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

if t<ConYawData(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(ConYawData(1,1)) ' instead.']);
    t = ConYawData(1,1);
elseif t>ConYawData(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(ConYawData(end,1)) ' instead.']);
    t = ConYawData(end,1);
end

yaw_o = interp1(ConYawData(:,1),ConYawData(:,2:end),t);
yaw = yaw_o(iT);
end

