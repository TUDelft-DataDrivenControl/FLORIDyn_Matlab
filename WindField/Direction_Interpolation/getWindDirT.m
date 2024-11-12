function phi = getWindDirT(WindDir,iT,t)
%GETWINDDIRT returns the wind direction at the respective turbine(s)
% Uniform interpolation version - all turbines experience the same changes
%   Expects a .csv called "WindTI" with the structure
%       time phi
%       time phi
%        ...
%       time phi
%   The value is interpolated linearly between the setpoints
% ======================================================================= %
% iT        = single value or array with turbine index / indices
% WindVel   = (t,phi) pairs between which is linearly interpolated 
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

% SOWFA angle
phi = ones(size(iT))*interp1(WindDir(:,1),WindDir(:,2),t);
end

