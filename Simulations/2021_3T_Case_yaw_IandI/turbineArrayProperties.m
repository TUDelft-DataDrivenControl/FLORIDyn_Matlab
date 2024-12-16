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

% Foot position of the turbine x,y,z (NOT nacelle position, Ground level!)
tProp.Pos = [...
            608  500 0;...     % T0
            1500 500 0;...     % T1
            2392 500 0;...     % T2
            ]; 

% Turbine type, the parameters are stored centralized in
% Data/TurbineData/getTurbineData.m
tProp.Type = {...
            'DTU 10MW',...
            'DTU 10MW',...
            'DTU 10MW',...
            };
      
% Initial world orientation of the Wind Turbine. Direction is 
% Orthogonal for wind from

% Has to fit the turbine states structure of the selected FLORIS model (see
% States.m), Gaussian: 'a','yaw'(deg,counter-clockwise,relative),'TI'
tProp.Init_States = [...
    0.33 0 0.06; ...
    0.33 0 0.06; ...
    0.33 0 0.06;
    ];

% NOTE ON THE TURBINE ORIENTATION, THE YAW ANGLE AND THE WIND DIRECTION
%   The wind direction and the turbine orientation are clockwise defined,
%   the yaw angle is the result of 
%       wind direction - orientation = yaw angle
%   Thus, the yaw angle is defined in clockwise direction.
%   Wind from -x to +x has the wind direction 270 deg.
%   Wind from +y to -y has the wind direction 0 / 360 deg.
end