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

function WSE = WSEParameters(nT,pathSimulation,deltaT)
%% Settings for the I & I wind speed estimator
try
    rotorSpeed  = ...
        importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_rotorSpeed.csv']);
catch
    rotorSpeed  = ...
        importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_rotorSpeedFiltered.csv']);
end
bladePitch  = importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_bladePitch.csv']);
genTorque   = importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_generatorTorque.csv']);
nacelleYaw  = importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_nacelleYaw.csv']);

% Set parameters
WSE.V       = ones(nT,1)*8; % Initial wind speed
WSE.Vinit   = WSE.V;
WSE.gamma   = 5;
WSE.beta    = 0;
WSE.omega   = rotorSpeed((1:nT),3)* pi/30;
WSE.Ee      = ones(nT,1)*0;

% Time during which the Initial wind speed will be used while the Estimator
% converges
WSE.Offset  = 100; % [s]

% Save turbine and Sim properties
WSE.T_prop  = estimator_dtu10mw();
WSE.dt_SOWF = rotorSpeed(nT+1,2)-rotorSpeed(1,2);
WSE.dt_FDyn = deltaT;
WSE.nT      = nT;

% Save data
WSE.rotorSpeed = rotorSpeed;
WSE.bladePitch = bladePitch;
WSE.genTorque  = genTorque;
WSE.nacelleYaw = nacelleYaw;

end