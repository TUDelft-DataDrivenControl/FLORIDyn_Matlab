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
function tProp = turbineArrayProperties()
%% Return location and type of the turbines

% Can be rotated around 2500, 2500 in a 5x5km domain with always
% maintaining at least 2.8D to an inflow plane.
tProp.Pos = [
1347.6   919.0 0.0
1640.9	1662.5 0.0
2248.0	1001.0 0.0
1934.2	2406.0 0.0
3149.3	1083.9 0.0
2227.5	3149.5 0.0
4218.2	1564.4 0.0
2520.8	3893.0 0.0
4097.0	2666.3 0.0
4036.4	3711.2 0.0];

tProp.Type = {...
            'DTU 10MW',...
            'DTU 10MW',...
            'DTU 10MW',...
            'DTU 10MW',...
            'DTU 10MW',...
            'DTU 10MW',...
            'DTU 10MW',...
            'DTU 10MW',...
            'DTU 10MW',...
            'DTU 10MW'};

% Initial world orientation of the Wind Turbine. Direction is 
% Orthogonal for wind from

% Has to fit the turbine states structure of the selected FLORIS model (see
% States.m), Gaussian: 'a','yaw'(deg,counter-clockwise,relative),'TI'
tProp.Init_States = [...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    ];

% NOTE ON THE TURBINE ORIENTATION, THE YAW ANGLE AND THE WIND DIRECTION
%   The wind direction and the turbine orientation are clockwise defined,
%   the yaw angle is the result of 
%       wind direction - orientation = yaw angle
%   Thus, the yaw angle is defined in clockwise direction.
%   Wind from -x to +x has the wind direction 270 deg.
%   Wind from +y to -y has the wind direction 0 / 360 deg.
end