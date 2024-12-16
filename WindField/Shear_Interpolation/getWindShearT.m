
% Copyright (C) <2024>, M Becker
%
% List of the contributors to the development of FLORIDyn: see LICENSE file.
% Description and complete License: see LICENSE file.
	
% This program (FLORIDyn) is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program (see COPYING file).  If not, see <https://www.gnu.org/licenses/>.
% ======================================================================= %
% Updated: 16. Dez. 2024, M. Becker
% ======================================================================= %

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

