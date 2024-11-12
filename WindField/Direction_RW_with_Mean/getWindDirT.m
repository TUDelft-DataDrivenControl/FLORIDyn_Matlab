function phi = getWindDirT(WindDirNow,WindDir)
%GETWINDDIRT returns the wind direction at the respective turbine(s)
% ======================================================================= %
% Random walk model with mean implementation
% ======= Input ======
% WindDir.Init   = Mean direction
% WindDir.ColSig = Tn x Tn = chol(Covariance Matrix)
% WindDirNow     = Current value
% iT             = Index/Indeces of the turbines
% t              = time of request
% ======================================================================= %
weigthedRandN = randn(1,length(WindDirNow));
phi = WindDirNow + (weigthedRandN*WindDir.CholSig)' + ...
    WindDir.MeanPull*(WindDir.Init - WindDirNow);
end