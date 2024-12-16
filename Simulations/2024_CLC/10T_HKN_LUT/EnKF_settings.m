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
function EnKF = EnKF_settings()
%ENKF_SETTINGS Returns a struct with the relevant Ensemble Kalman Filter
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
%   was found to be optimal for the EnKF by Whitaker and Hamill (2002).
% ======================================================================= %
%% Parameter tuning
% Best sample
    %   FLORIDyn    \sig u,dw 0.6966  \sig u,cw  0.3570 \sig u,t 206.2331    0.7396    
    %   FLORIS      pp 2.5433  alpha  0.7515  beta  0.3483  ka  0.7701  kfb  0.0467  kfc  0.2055  kTI  1.5374    
    %   EnKF        loc u 6.8011*D    pn_u 0.1991    mn_p 0.2633 (prev. 1000, 0.4, 0.1)
%% Settings
% Store results
EnKF.StoreData = false;
%% Simulation
% Number of Ensembles
EnKF.nE = 50;

% Number of steps between the corrections. (Integer)
%   A lower number will lead to a higher computational load, a larger 
%   number will lead to more diverging ensembles.
%   e.g. EnKF.nS = 4; will lead to a correction every fourth step.
%       NOTE: If the number of sim steps is not divisible by nS, the
%       simulation will be shortend to a fitting number of steps.
EnKF.nS = 3;

%% ================ Correction
%% Velocity
EnKF.Vel.Correct        = true;
%   Locatization
EnKF.Vel.loc            = true;
EnKF.Vel.cutOffLength   = 6.8011 * 178.4; %m
%   Inflation
EnKF.Vel.inf            = false;
EnKF.Vel.infFactor      = 1.07; 
%   Prediction model of how the wind field state will change during the
%   simulation. 'RW' is a random walk model, 'ZOH' keeps the most recent 
%   value. 'Interpolate' interpolates between given setpoints.
EnKF.Vel.PredicionModel = 'RW'; 
%   The state is polluted by process noise which mirrors the potential
%   natural development of the state. Is applied at every time step of the
%   simulation (potentially higher feq. than EnKF corrections!)
EnKF.Vel.ProcessNoise   = 0.1991; %0.4;%1.5;%0.4;% m/s 0.12; % m/s0.5;%

%   Measurement noise of the power generated
EnKF.Output.PowNoiseVar = 0.08;%0.2633;%0.1;%0.1%0.05;  % MW

%   Limits
%       Corrects all states to be within the min/max range - useful to
%       limit values to physical domain.
EnKF.Vel.Limits         = true;
EnKF.Vel.Min            = 1; %m/s
EnKF.Vel.Max            = 20; %m/s


%% Direction
EnKF.Dir.Correct        = true;
%   Locatization
EnKF.Dir.loc            = true;
EnKF.Dir.cutOffLength   = 500; %m
%   Inflation
EnKF.Dir.inf            = false;
EnKF.Dir.infFactor      = 1.07; 
%   Prediction model of how the wind field state will change during the
%   simulation. 'RW' is a random walk model, 'ZOH' keeps the most recent 
%   value. 'Interpolate' interpolates between given setpoints.
EnKF.Dir.PredicionModel = 'RW'; 
%   The state is polluted by process noise which mirrors the potential
%   natural development of the state. Is applied at every time step of the
%   simulation (potentially higher feq. than EnKF corrections!)
EnKF.Dir.ProcessNoise   = 3; % deg %1

EnKF.Dir.MeasurementNoise = 3;%.2; % Deg
%% Turbulence Intensity
EnKF.TI.Correct         = false;
%   Locatization
EnKF.TI.loc             = false;
EnKF.TI.cutOffLength    = 1000; %m
% Inflation
EnKF.TI.inf             = false;
EnKF.TI.infFactor       = 1.07; 
%   Prediction model of how the wind field state will change during the
%   simulation. 'RW' is a random walk model, 'ZOH' keeps the most recent 
%   value. 'Interpolate' interpolates between given setpoints.
EnKF.TI.PredicionModel  = 'ZOH'; 
%   The state is polluted by process noise which mirrors the potential
%   natural development of the state. Is applied at every time step of the
%   simulation (potentially higher feq. than EnKF corrections!)
EnKF.TI.ProcessNoise    = 0.005; % (percent)

%% Visulalization
% Set a fixed color palette to ensure that the steps are coloured equally.
EnKF.Colors = [...
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
EnKF.nColors = size(EnKF.Colors,1);

end