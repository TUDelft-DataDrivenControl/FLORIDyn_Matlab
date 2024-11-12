function phi = getWindDirT(WindDir,iT,t)
%GETWINDDIRT Return wind direction in SOWFA-deg for the requested
%turbine(s)
% ======================================================================= %
% Individual Turbine value implementation
%   requires a .csv in the simulation folder called WindDirTurbine.csv
%   where each row is a
%       time, phi_T0, phi_T1, ... phi_Tn
%   setpoint in time. The values are interploated linearly between the
%   setpoints.
% ======= Input ======
% WindDir   = (t,phi_T0, phi_T1, ... phi_Tn)
% iT        = Index/Indeces of the turbines
% t         = time of request
% ======================================================================= %
if t<WindDir(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindDir(1,1)) ' instead.']);
    t = WindDir(1,1);
elseif t>WindDir(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindDir(end,1)) ' instead.']);
    t = WindDir(end,1);
end

phi_out = interp1(WindDir(:,1),WindDir(:,2:end),t);
phi = phi_out(iT)';
end

