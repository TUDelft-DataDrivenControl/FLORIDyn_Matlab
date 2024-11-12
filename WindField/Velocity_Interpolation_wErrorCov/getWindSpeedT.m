function Vel = getWindSpeedT(WindVel,iT,t)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% Uniform interpolation version - all turbines experience the same changes
%   Expects a .csv called "WindTI" with the structure
%       time u
%       time u
%        ...
%       time u
%   The value is interpolated linearly between the setpoints
% ======================================================================= %
% WindVel.Data   = (t,u) pairs between which is linearly interpolated 
% WindVel.ColSig = Tn x Tn = col(Covariance Matrix)
% iT             = single value or array with turbine index / indices
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

% Set data and pollute with a random given covariance error
Vel = ones(size(iT))*interp1(WindVel.Data(:,1),WindVel.Data(:,2),t);
Vel = Vel + (randn(1,length(Vel))*WindVel.CholSig)';
end

