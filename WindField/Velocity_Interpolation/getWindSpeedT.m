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

function U = getWindSpeedT(WindVel,iT,t)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% Uniform interpolation version - all turbines experience the same changes
%   Expects a .csv called "WindTI" with the structure
%       time u
%       time u
%        ...
%       time u
%   The value is interpolated linearly between the setpoints
% ======================================================================= %
% iT        = single value or array with turbine index / indices
% WindVel   = (t,u) pairs between which is linearly interpolated 
% t         = time of request
% ======================================================================= %
if t<WindVel(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindVel(1,1)) ' instead.']);
    t = WindVel(1,1);
elseif t>WindVel(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindVel(end,1)) ' instead.']);
    t = WindVel(end,1);
end

U = ones(size(iT))*interp1(WindVel(:,1),WindVel(:,2),t);
end

