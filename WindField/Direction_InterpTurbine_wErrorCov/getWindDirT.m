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
% WindDir.Data   = (t,phi_T0, phi_T1, ... phi_Tn)
% WindDir.ColSig = nT x nT, col(Covariance Matrix)
% iT        = Index/Indeces of the turbines
% t         = time of request
% ======================================================================= %
if t<WindDir.Data(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindDir.Data(1,1)) ' instead.']);
    t = WindDir.Data(1,1);
elseif t>WindDir.Data(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindDir.Data(end,1)) ' instead.']);
    t = WindDir.Data(end,1);
end

phi_out = interp1(WindDir.Data(:,1),WindDir.Data(:,2:end),t);
phi = phi_out(iT);
phi = phi + (randn(1,length(phi))*WindDir.CholSig)';
end

