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

function Vel = getWindSpeedT(WindVel,iT,t)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% ======================================================================= %
% Individual Turbine value implementation
%   requires a .csv in the simulation folder called WindVelTurbine.csv
%   where each row is a
%       time, U_T0, U_T1, ... U_Tn
%   setpoint in time. The values are interploated linearly between the
%   setpoints.
% ======= Input ======
% WindVel.Data   = (t,U_T0, U_T1, ... U_Tn)
% WindVel.ColSig = Tn x Tn = col(Covariance Matrix)
% iT             = Index/Indeces of the turbines
% t              = time of request
% ======================================================================= %
if t<WindVel.Data(1,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindVel.Data(1,1)) ' instead.']);
    t = WindVel.Data(1,1);
elseif t>WindVel.Data(end,1)
    warning(['The time ' num2str(t) ' is out of bounds, will use '...
        num2str(WindVel.Data(end,1)) ' instead.']);
    t = WindVel.Data(end,1);
end

WindVel_out = interp1(WindVel.Data(:,1),WindVel.Data(:,2:end),t);

% Extract needed velocity and pollute with a random given covariance error
Vel = WindVel_out(iT);
Vel = Vel + (randn(1,length(Vel))*WindVel.CholSig)';
end

