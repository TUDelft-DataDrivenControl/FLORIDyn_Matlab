function U = getWindSpeedT(WindVel,iT,~)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% iT        = single value or array with turbine index / indices
% WindVel   = Depends on the wind speed method chosen

U = WindVel*ones(size(iT));
end

