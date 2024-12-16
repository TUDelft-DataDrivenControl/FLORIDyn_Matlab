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

%% Add paths
% Reset paths
restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED

%disp(['Current folder:' pwd])
% Add basic paths
addpath(genpath(['.' filesep 'Visualization']));
addpath(['.' filesep 'FLORIDynCL']);
addpath(['.' filesep 'Data' filesep 'TurbineData']);
addpath(['.' filesep 'Data' filesep 'StoreData']);
addpath(['.' filesep 'Correction' filesep 'Functions'])
addpath(['.' filesep 'EnsembleKalmanFilter'])
addpath(['.' filesep 'ClosedLoopControl'])
addpath(['.' filesep 'Correction' filesep 'GetData'])

% Load simulation folder
addpath(genpath(['.' filesep 'Simulations' filesep pathToSimulation]));
rmpath(genpath(['.' filesep 'Simulations' filesep pathToSimulation filesep 'Results']))