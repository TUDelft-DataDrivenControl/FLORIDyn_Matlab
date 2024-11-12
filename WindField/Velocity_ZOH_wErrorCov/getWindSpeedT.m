function Vel = getWindSpeedT(Vel, WindVelCholSig)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% WindVel        = previous wind speed 
% WindVelColSig  = Tn x Tn = col(Covariance Matrix)
% iT             = single value or array with turbine index / indices

Vel = Vel + (randn(1,length(Vel))*WindVelCholSig)';
end

