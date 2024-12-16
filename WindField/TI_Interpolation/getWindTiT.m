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

function Ti = getWindTiT(WindTi,iT,t)
%GETWINDTI Return turbulence intensity for the requested turbine(s)
% Uniform interpolation version - all turbines experience the same changes
%   Expects a .csv called "WindTI" with the structure
%       time TI
%       time TI
%        ...
%       time TI
%   The value is interpolated linearly between the setpoints
% ======= Input ======
% WindTi    = (t,TI) pairs between which is linearly interpolated
% iT        = Index/Indeces of the turbines
% t         = time of request
% ======================================================================= %

if t<WindTi(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindTi(1,1)) ' instead.']);
    t = WindTi(1,1);
elseif t>WindTi(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindTi(end,1)) ' instead.']);
    t = WindTi(end,1);
end

Ti = ones(size(iT))*interp1(WindTi(:,1),WindTi(:,2),t);
end

