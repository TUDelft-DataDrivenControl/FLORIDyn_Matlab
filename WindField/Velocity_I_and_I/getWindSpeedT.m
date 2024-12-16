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

function [U,WindVel] = getWindSpeedT(WindVel,iT,SimTime,WindDir,p_p)
%GETWINDSPEEDT returns the wind speed at the respective turbine(s)
% iT        = single value or array with turbine index / indices
% Wind      = Data needed for the I&I wind speed estimator
% SimTime   = Current simulation time
% =====================================================
% U         = Velocity
% WindVel   = Data needed for the I&I wind speed estimator

if SimTime == WindVel.StartTime
    U = WindVel.WSE.V(iT);
    return
end

% Calculate index of current time step and previous one
indCurr = round((SimTime-WindVel.StartTime)/WindVel.WSE.dt_SOWF);
indPrev = (WindVel.TimePrev-WindVel.StartTime)/WindVel.WSE.dt_SOWF+1;

% Number of Turbines
nT      = WindVel.WSE.nT;

% Run I&I from last time step to current one
for i = indPrev:indCurr
    % Get input data from SOWFA
    try
        Rotor_Speed = WindVel.WSE.rotorSpeed((1:nT) + nT*(i-1),3);
    catch
        dips('Failed to get Wind Speed');
    end
    Blade_pitch = WindVel.WSE.bladePitch((1:nT) + nT*(i-1),3);
    Gen_Torque  = WindVel.WSE.genTorque((1:nT) + nT*(i-1),3);
    
    yaw = deg2rad(WindDir - WindVel.WSE.nacelleYaw((1:nT) + nT*(i-1),3));
    
    % Run estimator
    [V_out,WindVel.WSE] = ...
        WindSpeedEstimatorIandI_FLORIDyn(...
        WindVel.WSE, Rotor_Speed, Blade_pitch, Gen_Torque,yaw,p_p);
end

% Is equal to SimTime if SimTime is divisable by SOWFA time steps
%   If not, this way no steps are skipped
WindVel.TimePrev = indCurr*WindVel.WSE.dt_SOWF + WindVel.StartTime;
if (SimTime-WindVel.StartTime)>WindVel.WSE.Offset
    U = WindVel.WSE.V(iT);
else
    U = WindVel.WSE.Vinit(iT);
end
end

