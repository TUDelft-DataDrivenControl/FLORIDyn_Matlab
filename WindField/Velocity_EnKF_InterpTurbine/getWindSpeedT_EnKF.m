function Vel = getWindSpeedT_EnKF(WindVel,iT,t)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% ======================================================================= %
% Individual Turbine value implementation
%   requires a .csv in the simulation folder called WindVelTurbine.csv
%   where each row is a
%       time, U_T0, U_T1, ... U_Tn
%   setpoint in time. The values are interploated linearly between the
%   setpoints.
% ======= Input ======
% WindVel   = (t,U_T0, U_T1, ... U_Tn)
% iT        = Index/Indeces of the turbines
% t         = time of request
% ======================================================================= %
if t<WindVel(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindVel(1,1)) ' instead.']);
    t = WindVel(1,1);
elseif t>WindVel(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindVel(end,1)) ' instead.']);
    t = WindVel(end,1);
end

WindVel_out = interp1(WindVel(:,1),WindVel(:,2:end),t);
Vel = WindVel_out(iT);
end

