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

function UKF = UKF_settings()
%UKF_SETTINGS Returns a struct with the relevant Unscented Kalman Filter
%settings.
% ======================================================================= %
% Covariance Localization
%   If locatlization is used the influence of the measurements is
%   weighted by the distance to the location of the measurement. In this
%   simulation turbines are assumed to be the location of the measurement.
%   It is encouraged to use the feature.
%       (Houtekamer and Mitchell, 2001; Hamill et al., 2001; Whitaker
%        and Hamill, 2002)
%   The cut off length determines how far the correction reaches. This is
%   not equal to an influence of 0. No influence is reached at a distance
%   of 2 * sqrt(10/3) * cutOffLength. (Gaspari and Cohn (1999))
% ======================================================================= %
% Covariance inflation
%   This means that the state error is artificially inflated. It can help
%   to find a better solution, even though the solution is undersampled.
%   This in return allows the use of less ensembles.
%   Recommended values are 1,3 or 7% depending on the implementation. 7%
%   was found to be optimal for the UKF by Whitaker and Hamill (2002).
% ======================================================================= %
%% Settings
%% Simulation
% Number of sigma points (has to be uneven)
UKF.nE = 5;
UKF.L       = (UKF.nE-1)/2;
UKF.alpha   = 1e-3;         % from [1]
UKF.beta    = 2;            % from [1]
UKF.kappa   = 0;            % from [1]
UKF.lambda  = UKF.alpha^2*(UKF.L + UKF.kappa) - UKF.L; % from [1]
UKF.W0      = 0;            % from [2]

% Number of steps between the corrections. (Integer)
%   A lower number will lead to a higher computational load, a larger 
%   number will lead to more diverging ensembles.
%   e.g. UKF.nS = 4; will lead to a correction every fourth step.
%       NOTE: If the number of sim steps is not divisible by nS, the
%       simulation will be shortend to a fitting number of steps.
UKF.nS = 3;

%% ================ Correction
%% Velocity
UKF.Vel.Correct        = true;
%   Locatization
UKF.Vel.loc            = true;
UKF.Vel.cutOffLength   = 500; %m
%   Inflation
UKF.Vel.inf            = false;
UKF.Vel.infFactor      = 1.07; 
%   Prediction model of how the wind field state will change during the
%   simulation. 'RW' is a random walk model, 'ZOH' keeps the most recent 
%   value. 'Interpolate' interpolates between given setpoints.
UKF.Vel.PredicionModel = 'RW'; 
%   The state is polluted by process noise which mirrors the potential
%   natural development of the state. Is applied at every time step of the
%   simulation (potentially higher feq. than UKF corrections!)
UKF.Vel.ProcessNoise   = 0.1;%0.4;% m/s 0.12; % m/s

%   Measurement noise of the power generated
UKF.Output.PowNoiseVar = 0.05;%0.1%0.05;  % MW
%% Direction
UKF.Dir.Correct        = true;
%   Locatization
UKF.Dir.loc            = true;
UKF.Dir.cutOffLength   = 1000; %m
%   Inflation
UKF.Dir.inf            = false;
UKF.Dir.infFactor      = 1.07; 
%   Prediction model of how the wind field state will change during the
%   simulation. 'RW' is a random walk model, 'ZOH' keeps the most recent 
%   value. 'Interpolate' interpolates between given setpoints.
UKF.Dir.PredicionModel = 'RW'; 
%   The state is polluted by process noise which mirrors the potential
%   natural development of the state. Is applied at every time step of the
%   simulation (potentially higher feq. than UKF corrections!)
UKF.Dir.ProcessNoise   = 3; % deg %1

%% Turbulence Intensity
UKF.TI.Correct         = false;
%   Locatization
UKF.TI.loc             = false;
UKF.TI.cutOffLength    = 1000; %m
% Inflation
UKF.TI.inf             = false;
UKF.TI.infFactor       = 1.07; 
%   Prediction model of how the wind field state will change during the
%   simulation. 'RW' is a random walk model, 'ZOH' keeps the most recent 
%   value. 'Interpolate' interpolates between given setpoints.
UKF.TI.PredicionModel  = 'ZOH'; 
%   The state is polluted by process noise which mirrors the potential
%   natural development of the state. Is applied at every time step of the
%   simulation (potentially higher feq. than UKF corrections!)
UKF.TI.ProcessNoise    = 0.005; % (percent)

%% Visulalization
% Set a fixed color palette to ensure that the steps are coloured equally.
UKF.Colors = [...
     12,  35,  64 
    255, 184,  28 
    108, 194,  74
      0, 118, 194
    237, 104,  66 
    0,   155, 119 
    224,  60,  49
    111,  29, 119 
    0,   184, 200 
    239,  96, 163 
	165,   0,  52]./255;
UKF.ColIFAC = [...
     96 142 185     % Blue
    230   0  19     % Red
    191 157  90     % Gold
    219 214 214     % Warm white
    ]./255;
UKF.nColors = size(UKF.Colors,1);

end
% [1] The Unscented Kalman Filter for Nonlinear Estimation
%   Eric A. Wan and Rudolph van der Merwe
% [2] Unscented Kalman Filter Tutorial
%   Gabriel A. Terejanu