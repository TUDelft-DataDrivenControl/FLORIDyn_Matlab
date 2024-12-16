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

function states = States()
%STATES States of the FLORIS model
%   T  -> turbine states needed to calculate the wake
%   OP -> Observation point states (world and wake coordinates)
%   WF -> Wind field states
% ======================================================================= %
states.T_names  = {'a','yaw','TI'};
states.Turbine  = length(states.T_names);

states.OP_names = {'x0','y0','z0','x1','y1','z1'};
states.OP       = length(states.OP_names);

states.WF_names = {'wind_vel','wind_dir','TI0'};
states.WF       = length(states.WF_names);
end

