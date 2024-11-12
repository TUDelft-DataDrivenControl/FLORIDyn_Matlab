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
% WindDir.Data   = (t,phi) pairs between which is linearly interpolated 
% WindDir.ColSig = nT x nT, col(Covariance Matrix)
% iT        = single value or array with turbine index / indices
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

% SOWFA angle
phi = ones(size(iT))*interp1(WindDir.Data(:,1),WindDir.Data(:,2),t);
phi = phi + (randn(1,length(phi))*WindDir.CholSig)';
end

