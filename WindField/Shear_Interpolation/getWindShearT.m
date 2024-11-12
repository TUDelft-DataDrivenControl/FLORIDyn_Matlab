function shear = getWindShearT(WindShear,z)
%GETWINDSHEART returns the relative reduction (or speed-up) by wind shear
% Expects a .csv file called "WindShearProfile.csv" with a normalized wind
% speed profile for different heights:
%   z, (u_z/u0)
%   z, (u_z/u0)
%   z, (u_z/u0)
% There is a linear interpolation between every pair
% IN CASE z IS OUT OF BOUNDS the function will use the closest available
% setpoint
% ======================================================================= %
% WindShear = normalized wind speed at different heights
% z         = height(s)
% ======================================================================= %
% Out of bounds handling
maxZ = max(WindShear(:,1));
minZ = min(WindShear(:,1));
z(z>maxZ) = maxZ;
z(z<minZ) = minZ;
% Interpolate
shear = interp1(WindShear(:,1),WindShear(:,2),z);
end

