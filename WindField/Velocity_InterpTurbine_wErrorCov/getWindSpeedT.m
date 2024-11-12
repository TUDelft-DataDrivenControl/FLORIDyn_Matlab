function Vel = getWindSpeedT(WindVel,iT,t)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% ======================================================================= %
% Individual Turbine value implementation
%   requires a .csv in the simulation folder called WindVelTurbine.csv
%   where each row is a
%       time, U_T0, U_T1, ... U_Tn
%   setpoint in time. The values are interploated linearly between the
%   setpoints.
% ======= Input ======
% WindVel.Data   = (t,U_T0, U_T1, ... U_Tn)
% WindVel.ColSig = Tn x Tn = col(Covariance Matrix)
% iT             = Index/Indeces of the turbines
% t              = time of request
% ======================================================================= %
if t<WindVel.Data(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindVel.Data(1,1)) ' instead.']);
    t = WindVel.Data(1,1);
elseif t>WindVel.Data(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindVel.Data(end,1)) ' instead.']);
    t = WindVel.Data(end,1);
end

WindVel_out = interp1(WindVel.Data(:,1),WindVel.Data(:,2:end),t);

% Extract needed velocity and pollute with a random given covariance error
Vel = WindVel_out(iT);
Vel = Vel + (randn(1,length(Vel))*WindVel.CholSig)';
end

