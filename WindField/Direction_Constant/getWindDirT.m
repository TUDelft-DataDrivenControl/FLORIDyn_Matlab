function phi = getWindDirT(WindDir,iT,~)
%GETWINDDIRT Return wind direction in SOWFA-deg for the requested
%turbine(s)
%
% ======= Input ======
% WindDir   = 
% iT        = Index/Indeces of the turbines

phi = ones(size(iT))*WindDir;
end

