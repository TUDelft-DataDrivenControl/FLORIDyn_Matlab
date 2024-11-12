function yaw = getYaw(ConYawData,iT,t)
%GETYAW Return the yaw angle for the requested turbine(s)
% ======================================================================= %
% SOWFA implementation
%   requires a .csv in the simulation folder called SOWFA_nacelleYaw.csv
%   from a SOWFA simulation
% ======= Input ======
% ConYawData = SOWFA data read by importSOWFAFile()
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

