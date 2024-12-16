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

%% Links to the specific version of FLORIS which is used
% NOTE: There are no values to be modified in this file, it has to be in
% the simulation folder to automatically retreive the path to the
% simulation folder and store it in "Sim.PathToSim".
%
% ======= DO NOT MODIFY ================================================= %
if ispc
    % Microsoft pathing
    addpath(['.\FLORIS\' Sim.FLORIS]);
    addpath(['.\WindField\Velocity_' Wind.Input_Vel]);
    addpath(['.\WindField\Direction_' Wind.Input_Dir]);
    addpath(['.\WindField\TI_' Wind.Input_TI]);
    addpath(['.\WindField\Shear_' Wind.Input_Shr]);
    addpath(['.\FLORIS\Discretization\' Sim.RotorDiscret]);
    addpath(['.\Controller\Yaw_' Con.Yaw]);
    addpath(['.\Correction\Direction_' Wind.Correction.Dir]);
    addpath(['.\Correction\Velocity_' Wind.Correction.Vel]);
    addpath(['.\Correction\TI_' Wind.Correction.TI]);
else
    % Unix pathing
    addpath(['FLORIS/' Sim.FLORIS]);
    addpath(['WindField/Velocity_' Wind.Input_Vel]);
    addpath(['WindField/Direction_' Wind.Input_Dir]);
    addpath(['WindField/TI_' Wind.Input_TI]);
    addpath(['WindField/Shear_' Wind.Input_Shr]);
    addpath(['FLORIS/Discretization/' Sim.RotorDiscret]);
    addpath(['Controller/Yaw_' Con.Yaw]);
    addpath(['Correction/Direction_' Wind.Correction.Dir]);
    addpath(['Correction/Velocity_' Wind.Correction.Vel]);
    addpath(['Correction/TI_' Wind.Correction.TI]);
end
p = mfilename('fullpath');
Sim.PathToSim = p(1:end-14);
clear p
% ======= DO NOT MODIFY ================================================= %