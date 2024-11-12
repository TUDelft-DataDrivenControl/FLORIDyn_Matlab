function Vel = getWindSpeedT(WindVel,iT,~)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% WindVel.Data   = u value 
% WindVel.ColSig = Tn x Tn = col(Covariance Matrix)
% iT             = single value or array with turbine index / indices

Vel = WindVel.Data*ones(size(iT));
Vel = Vel + (randn(1,length(Vel))*WindVel.CholSig)';
end

