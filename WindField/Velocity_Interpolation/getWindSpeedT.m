function U = getWindSpeedT(WindVel,iT,t)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% Uniform interpolation version - all turbines experience the same changes
%   Expects a .csv called "WindTI" with the structure
%       time u
%       time u
%        ...
%       time u
%   The value is interpolated linearly between the setpoints
% ======================================================================= %
% iT        = single value or array with turbine index / indices
% WindVel   = (t,u) pairs between which is linearly interpolated 
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

U = ones(size(iT))*interp1(WindVel(:,1),WindVel(:,2),t);
end

