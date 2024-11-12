function phi = getWindDirT(WindDir,iT,~)
%GETWINDDIRT Return wind direction in SOWFA-deg for the requested
%turbine(s)
%
% ======= Input ======
% WindDir.Data   = float, wind direction
% WindDir.ColSig = nT x nT, col(Covariance Matrix)
% iT        = Index/Indeces of the turbines

phi = ones(size(iT))*WindDir.Data;
phi = phi + (randn(1,length(phi))*WindDir.CholSig)';
end

