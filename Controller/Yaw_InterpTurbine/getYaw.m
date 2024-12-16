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

function yaw = getYaw(ConYawData,iT,t)
%GETYAW Return the yaw angle for the requested turbine(s)
% ======================================================================= %
% Individual Turbine value implementation
%   requires a .csv in the simulation folder called Control_YawInterp.csv
%   where each row is a
%       time, yaw_T0, yaw_T1, ... yaw_Tn
%   setpoint in time. The values are interploated linearly between the
%   setpoints.
% ======= Input ======
% ConYawData = (t,yaw_T0, yaw_T1, ... yaw_Tn)
% iT         = Index/Indeces of the turbines
% t          = time of request
% ======================================================================= %

if t<ConYawData(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(ConYawData(1,1)) ' instead.']);
    t = ConYawData(1,1);
elseif t>ConYawData(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(ConYawData(end,1)) ' instead.']);
    t = ConYawData(end,1);
end

yaw_o = interp1(ConYawData(:,1),ConYawData(:,2:end),t);
yaw = yaw_o(iT);
end

