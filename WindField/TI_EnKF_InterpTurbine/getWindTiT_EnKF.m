function Ti = getWindTiT_EnKF(WindTi,iT,t)
%GETWINDTI Return turbulence intensity for the requested turbine(s)
% ======================================================================= %
% Individual Turbine value implementation
%   requires a .csv in the simulation folder called WindTITurbine.csv
%   where each row is a
%       time, TI_T0, TI_T1, ... TI_Tn
%   setpoint in time. The values are interploated linearly between the
%   setpoints.
% ======= Input ======
% WindTi    = (t,TI_T0, TI_T1, ... TI_Tn)
% iT        = Index/Indeces of the turbines
% t         = time of request
% ======================================================================= %

if t<WindTi(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindTi(1,1)) ' instead.']);
    t = WindTi(1,1);
elseif t>WindTi(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindTi(end,1)) ' instead.']);
    t = WindTi(end,1);
end
Ti_out = interp1(WindTi(:,1),WindTi(:,2:end),t);
Ti = Ti_out(iT);
end

