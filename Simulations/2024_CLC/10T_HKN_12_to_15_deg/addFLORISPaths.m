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
addpath(['.' filesep 'FLORIS' filesep Sim.FLORIS]);
addpath(['.' filesep 'WindField' filesep 'Velocity_' Wind.Input_Vel]);
addpath(['.' filesep 'WindField' filesep 'Direction_' Wind.Input_Dir]);
addpath(['.' filesep 'WindField' filesep 'TI_' Wind.Input_TI]);
addpath(['.' filesep 'WindField' filesep 'Shear_' Wind.Input_Shr]);
addpath(['.' filesep 'FLORIS' filesep 'Discretization' filesep Sim.RotorDiscret]);
addpath(['.' filesep 'Controller' filesep 'Yaw_' Con.Yaw]);
addpath(['.' filesep 'Correction' filesep 'Direction_' Wind.Correction.Dir]);
addpath(['.' filesep 'Correction' filesep 'Velocity_' Wind.Correction.Vel]);
addpath(['.' filesep 'Correction' filesep 'TI_' Wind.Correction.TI]);
p = mfilename('fullpath');
Sim.PathToSim = p(1:end-14);
clear p