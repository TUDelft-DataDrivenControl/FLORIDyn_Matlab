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

function [UKF, Wind] = UKF_PreProcessing(UKF, T, Wind, Sim)
%UKF_PREPROCESSING Summary of this function goes here
%   Detailed explanation goes here
%% Simulation Sections
UKF.Sim.Sections   = max(floor(Sim.nSimSteps / UKF.nS),1);
UKF.Sim.SecDur     = Sim.TimeStep * UKF.nS;
UKF.Sim.StartTime  = Sim.StartTime;

%% Ensemble States

% States
UKF.States.Vel = repmat(T.States_WF(:,1),1,UKF.nE);
UKF.States.Dir = repmat(T.States_WF(:,2),1,UKF.nE);

% Measurement noise (NOT process noise!)
%   Velocity
[UKF.Vel.C_ee_Vel, UKF.Vel.C_ee_Vel_Chol] = readCovMatrix(...
    readmatrix('WindVelCovariance.csv'),T.nT,'WindVel');
%   Direction
[UKF.Dir.C_ee_Dir, UKF.Dir.C_ee_Dir_Chol] = readCovMatrix(...
    readmatrix('WindDirCovariance.csv'),T.nT,'WindDir');
%   Power generated
[UKF.Output.C_ee_Pow, UKF.Output.C_ee_Pow_Chol] = readCovMatrix(...
    UKF.Output.PowNoiseVar,T.nT,'WindVel');

% ============= DISTRIBUTE SIGMA INITIAL POINTS ========================= %
% Eq.(15) in [1]
UKF.P = [UKF.Vel.ProcessNoise, 0; 0, UKF.Dir.ProcessNoise];
sp_off = sqrt((UKF.L + UKF.lambda)*UKF.P); 
for iE = 1:UKF.L
    UKF.States.Vel(:,iE+1) = UKF.States.Vel(:,iE+1) + sp_off(iE,1);
    UKF.States.Dir(:,iE+1) = UKF.States.Dir(:,iE+1) + sp_off(iE,2);
end

for iE = UKF.L+1:UKF.nE-1
    UKF.States.Vel(:,iE+1) = UKF.States.Vel(:,iE+1) - sp_off(iE-UKF.L,1);
    UKF.States.Dir(:,iE+1) = UKF.States.Dir(:,iE+1) - sp_off(iE-UKF.L,2);
end


% ============ Random approach ========================================== %
% % Velocity states
% UKF.States.Vel   = repmat(T.States_WF(:,1),1,UKF.nE);
% 
% % Velocity measurement noise (NOT process noise!)
% UKF.Vel.C_ee_Vel = readmatrix('WindVelCovariance.csv');
% [UKF.Vel.C_ee_Vel, UKF.Vel.C_ee_Vel_Chol] = readCovMatrix(...
%     readmatrix('WindVelCovariance.csv'),T.nT,'WindVel');
% for iE = 1:UKF.nE
%     randomError = (randn(T.nOP,T.nT)*EnKF.Vel.C_ee_Vel_Chol);
%     UKF.States.Vel(:,iE) = UKF.States.Vel(:,iE) + randomError (:);
% end
% 
% % Power generated noise
% [UKF.Output.C_ee_Pow, UKF.Output.C_ee_Pow_Chol] = readCovMatrix(...
%     UKF.Output.PowNoiseVar,T.nT,'WindVel');
% 
% 
% UKF.States.Dir = repmat(T.States_WF(:,2),1,UKF.nE);
% [UKF.Dir.C_ee_Dir, UKF.Dir.C_ee_Dir_Chol] = readCovMatrix(...
%     readmatrix('WindDirCovariance.csv'),T.nT,'WindDir');
% for iE = 1:UKF.nE
%     randomError = (randn(T.nOP,T.nT)*UKF.Dir.C_ee_Dir_Chol);
%     UKF.States.Dir(:,iE) = UKF.States.Dir(:,iE) + randomError (:);
% end


UKF.States_OP = repmat(T.States_OP,1,UKF.nE);
UKF.nStatesOP = length(T.Names_OP);

%% Measurements
UKF.M = cell(UKF.nE,1);
UKF.Interaction = cell(UKF.nE,T.nT);
UKF.InteractionNames = 'T_ID, iOP1, wOP1, iOP2, wOP2, w';
UKF.Output.Pow = zeros(T.nT,UKF.nE);

%% Flow field inputs
% Allocate fields
UKF.Wind = Wind;

%% Wind speed
switch lower(UKF.Vel.PredicionModel)
    case 'rw'
        Wind = rmfield(Wind,'Vel');
        Wind.Input_Vel = 'EnKF_RW';
        Wind.Vel.CholSig = UKF.Vel.C_ee_Vel_Chol;
    case 'zoh'
        Wind.Input_Vel = 'EnKF_ZOH';
        Wind.Vel = -1;
    case 'interpolate'
        Wind.Input_Vel = 'EnKF_InterpTurbine';
        Wind.Vel = [[...
            UKF.Sim.StartTime:Sim.TimeStep:UKF.Sim.StartTime + ...
            UKF.Sim.SecDur]',zeros(UKF.nS+1,T.nT)];
    otherwise
        Wind = rmfield(Wind,'Vel');
        Wind.Input_Vel = 'EnKF_RW';
        Wind.Vel.CholSig = UKF.Vel.C_ee_Vel_Chol;
end
% Enable process noise
%   Necessary to model the increasing uncertainty in areas where the
%   values can not be corrected and to encourange correction in areas
%   where values can be corrected.
Wind.Pertubation.Vel = true;
Wind.Pertubation.VelSigma = UKF.Vel.ProcessNoise;



%% Wind direction

switch lower(UKF.Dir.PredicionModel)
    case 'rw'
        Wind = rmfield(Wind,'Dir');
        Wind.Input_Dir = 'EnKF_RW';
        Wind.Dir.CholSig = UKF.Dir.C_ee_Dir_Chol;
    case 'zoh'
        Wind = rmfield(Wind,'Dir');
        Wind.Input_Dir = 'EnKF_ZOH';
        Wind.Dir = -1;
    case 'interpolate'
        Wind.Input_Dir = 'EnKF_InterpTurbine';
        Wind.Dir = [[...
            UKF.Sim.StartTime:Sim.TimeStep:UKF.Sim.StartTime + ...
            UKF.Sim.SecDur]',zeros(UKF.nS+1,T.nT)];
    otherwise
        Wind = rmfield(Wind,'Dir');
        Wind.Input_Dir = 'EnKF_RW';
        Wind.Dir.CholSig = UKF.Dir.C_ee_Dir_Chol;
end
% Enable process noise
%   Necessary to model the increasing uncertainty in areas where the
%   values can not be corrected and to encourange correction in areas
%   where values can be corrected.
Wind.Pertubation.Dir = true;
Wind.Pertubation.DirSigma = UKF.Dir.ProcessNoise;

end

% [1] The Unscented Kalman Filter for Nonlinear Estimation
%   Eric A. Wan and Rudolph van der Merwe