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

function Wind_Con = predict_wind(Wind, ~, ~, ~)
%PREDICT_WIND Takes the wind struct and generates a new one for the
%controller.
% ======================================================================= %
% INPUT (relevant)
%   Wind
%     .Input_Vel    (string)        Turbine measurement/input method
%     .Input_Dir    (string)        Turbine measurement/input method
%     .Input_TI     (string)        Turbine measurement/input method
%     .Correction
%       .Vel        (string)        OP state correction method
%       .Dir        (string)        OP state correction method
%       .TI         (string)        OP state correction method
%     .Vel          (misc. data)    Data for turbine measurements
%     .Dir          (misc. data)    Data for turbine measurements
%     .TI           (misc. data)    Data for turbine measurements
%     .Pertubation
%       .Vel        (boolean)       Pertubate OP states with random noise
%       .Dir        (boolean)       Pertubate OP states with random noise
%       .TI         (boolean)       Pertubate OP states with random noise
%       .VelSigma   (float)         Std of random noise
%       .DirSigma   (float)         Std of random noise
%       .TISigma    (float)         Std of random noise
%
%   Time            (float)         Current simulation time
%
%   CLC
%     .Time
%       .nS         (int)
%       .SecDur     (float)
%     .Con
%       .horizon_prediction (int)
%
%   deltaT          (float)         Seconds per time step

% ======================================================================= %
% Output
%   Wind_Con        Copy of Wind with adapted settings
%% ===================================================================== %%
% Keep wind speed and direction constant throught the simulation
Wind_Con = Wind;

%% Disable state pertubation
Wind_Con.Pertubation.Vel = false;
Wind_Con.Pertubation.Dir = false;
Wind_Con.Pertubation.TI  = false;

%% Set input mode
Wind_Con.Input_Vel = 'CLC_weighted_ZOH';
Wind_Con.Input_Dir = 'CLC_weighted_ZOH';
Wind_Con.Input_TI  = 'Constant';

%% Correction - Does not have an impact as the paths are already added
Wind_Con.Correction.Vel = 'None';
Wind_Con.Correction.Dir = 'None';
Wind_Con.Correction.TI  = 'None';

%% Generate data
% Not needed in this case
end

