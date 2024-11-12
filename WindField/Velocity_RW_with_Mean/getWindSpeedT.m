function Vel = getWindSpeedT(WindVelNow,WindVel)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% ======================================================================= %
% Random walk model with mean implementation
% ======= Input ======
% WindVel.Init    = Mean velocity
% WindVel.CholSig = Tn x Tn = col(Covariance Matrix)
% WindVelNow      = Current value
% iT              = Index/Indeces of the turbines
% t               = time of request
% ======================================================================= %
weigthedRandN = randn(1,length(WindVelNow));
Vel = WindVelNow + (weigthedRandN*WindVel.CholSig)' + ...
    WindVel.MeanPull*(WindVel.Init - WindVelNow);
end