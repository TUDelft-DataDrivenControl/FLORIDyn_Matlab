function Ti = getWindTiT(WindTi,iT,t)
%GETWINDTI Return turbulence intensity for the requested turbine(s)
% Uniform interpolation version - all turbines experience the same changes
%   Expects a .csv called "WindTI" with the structure
%       time TI
%       time TI
%        ...
%       time TI
%   The value is interpolated linearly between the setpoints
% ======= Input ======
% WindTi    = (t,TI) pairs between which is linearly interpolated
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

Ti = ones(size(iT))*interp1(WindTi(:,1),WindTi(:,2),t);
end

