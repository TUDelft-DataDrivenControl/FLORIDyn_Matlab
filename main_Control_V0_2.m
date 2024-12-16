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

%% FLORIDyn function for closed-loop-control with FLORIDyn
% Version 0.1
%   Features:
%       - Yaw control
%       - Uses FLORIDyn as optimization model
%       - Uses FLORIDyn as true model
%% ====== Choose simulation ======
pathToSimulation    = ['2024_CLC' filesep '10T_HKN_C1_201deg'];%'2023_Torque_BSA';%'2023_ACC_3T_yaw_steering';
pathToResults       = ['Results' filesep 'Control'];
caseName            = '2024_10T_HKN_CLC_001'; 
%% ====== Load data from the simulation ======
% Reset the Matlab Path and load essential paths & the simulation path
addPaths;

% Get the settings for the wind field, visualization, controller and Sim.
[Wind, Vis, Sim, Con] = setup();

% Add according functions to the search path
addFLORISPaths;
addFLORIDynPaths;
addEnKFPaths;

% Load linked data
turbProp        = turbineArrayProperties();
paramFLORIS     = parameterFLORIS();
paramFLORIDyn   = parameterFLORIDyn();
CLC             = clc_settings(Sim);

% _-_-_-_-_- OVERWRITE -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ %
% CLC.Con.Strategy            = Overwrite_Strategy;
% CLC.Con.DOFpT               = Overwrite_DOFpT;
% CLC.Con.horizon_prediction  = Overwrite_horizon_prediction;
% CLC.Con.horizon_action      = Overwrite_horizon_action;
% caseName                    = Overwrite_caseName;
% CLC.Time.Sections           = Overwrite_Sections;
% Sim.EndTime                 = Overwrite_EndTime;
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- %

addCLCPaths;
%%
turbProp.Init_States(:,2) = 0;
%% ====== Preprocess loaded data ======
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);
clear turbProp

%% ====== Init simulation ======
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

%% ====== Init closed loop parameters ======
CLC.gd = zeros(T.nT,1);
CLC.g0 = T.States_WF(T.StartI,2); % TODO Yaw_0 is WF & yaw

switch CLC.Con.InitMethod
    case "wind dir +- rand"
        CLC.x0 = mean(T.States_WF(T.StartI,2)) + ...
            (2*rand(CLC.Con.DOFpT*T.nT,1)-1);
    case "0 + rand"
        CLC.x0 = rand(CLC.Con.DOFpT*T.nT,1);
    otherwise
        CLC.x0 = rand(CLC.Con.DOFpT*T.nT,1);
end

Sim.nSimSteps = CLC.Time.nS;

Con.tanhYaw     = CLC.Con.tanhYaw;
Con.yawRangeMin = CLC.Con.yawRangeMin;
Con.yawRangeMax = CLC.Con.yawRangeMax;

% Power out
%   (:,1)           Time
%   (:,2)-(:,nT+1)  Turbine individual power
%   (:,nT+2)        Cumulative power
Pwr_out = zeros(CLC.Time.Sections*CLC.Time.nS,T.nT+2);

% Orientation applied
%   (:,1)           Time
%   (:,2)-(:,nT+1)  Turbine orientation
Ori_out = zeros(CLC.Time.Sections*CLC.Time.nS,T.nT+1);

% Orientation planned
%   (:,1)           Time
%   (:,2)-(:,nT+1)  Turbine orientation
Ori_pln = zeros(CLC.Time.Sections*(CLC.Con.horizon_action+1),T.nT+1); % USED TO BE + 2 DONT KNOW WHY

% Reduction of the free wind speed
%   (:,1)           Time
%   (:,2)-(:,nT+1)  percentages of the free wind speed perceived
Red_Uf = zeros(CLC.Time.Sections*CLC.Time.nS,T.nT+1);

% Interaction between the turbines
%   (:,1)           Time
%   (:,2)-(:,nT+1)  Reduction by the given turbine
Intrct = zeros(CLC.Time.Sections*CLC.Time.nS*T.nT,T.nT+1);

% Optimization parameters
%   (nS x nO)   Optimization value(s) over time
Opt_param = zeros(CLC.Time.Sections, CLC.Con.DOFpT*T.nT);

%% ====== Run simulation ======
for iCLC = 1:CLC.Time.Sections
    % \\\\ Core Loop \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % Estimation
    %   The from here on used state is considered the best estimate of the
    %   true system state
    
    % Optimize
    [Con, CLC] = controller(T, Wind, Sim, Con, Vis, paramFLORIDyn,...
        paramFLORIS,CLC, Sim.StartTime);
    
    % Run until the next control update
    [T, M, Vis, Mint] = ...
        FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
    
    
    % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


    % //// Store Variables ////////////////////////////////////////////////
    % Check if trajectories have been optimised, if not, skip storing vars
    if Sim.StartTime<CLC.Time.StartTime+0.001
        % Update time
        Sim.StartTime = Sim.StartTime + CLC.Time.SecDur;
        continue
    end

    %   Optimization parameters
    Opt_param(iCLC,:) = CLC.x0';

    %   Applied orientation
    Ori_out((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,:) = ...
        Con.YawData(1:CLC.Time.nS,:);

    %   Proposed orientation
    Ori_pln((1:(CLC.Con.horizon_action+1))+...          % USED TO BE + 2 DONT KNOW WHY
        (iCLC-1)*(CLC.Con.horizon_action+1),:) = ...    % USED TO BE + 2 DONT KNOW WHY
        Con.YawData(1:CLC.Con.horizon_action+1,:);      % USED TO BE + 2 DONT KNOW WHY
    
    %   Power & reduction timestamps
    Pwr_out((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,1) = ...
        M.("Time [s]")(1:T.nT:end);
    Red_Uf((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,1) = ...
        M.("Time [s]")(1:T.nT:end);

    %   Power turbines & reduction
    for iT = 1:T.nT
        Pwr_out((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,1+iT) = ...
            M.("Power generated [MW]")(iT:T.nT:end);
        Red_Uf((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,1+iT) = ...
            M.("Foreign Reduction [%]")(iT:T.nT:end);
    end
    
    % Interaction
    Intrct((1:CLC.Time.nS*T.nT)+(iCLC-1)*CLC.Time.nS*T.nT,:) = ...
        Mint;
    
    % /////////////////////////////////////////////////////////////////////
    Sim.StartTime = Sim.StartTime + CLC.Time.SecDur;

    
end
Pwr_out(:,end) = sum(Pwr_out(:,2:end-1),2);


%/\/\/\/\/\/\/\/\/\/\/\/\ Store outputs /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\%
folderName = [Sim.PathToSim pathToResults filesep caseName];
if not(isfolder(folderName))
    mkdir(folderName)
end
% Power out
storePath = [folderName filesep 'power_generated.csv'];
writematrix(Pwr_out,storePath)
% Orientation applied
storePath = [folderName filesep 'orientation_applied.csv'];
writematrix(Ori_out,storePath)
% Orientation planned
storePath = [folderName filesep 'orientation_planned.csv'];
writematrix(Ori_pln,storePath)
% Reduction of free wind speed
storePath = [folderName filesep 'wind_speed_reduction.csv'];
writematrix(Red_Uf,storePath)
% Interaction between the turbines
storePath = [folderName filesep 'interaction.csv'];
writematrix(Intrct,storePath)
% Optimization parameters
storePath = [folderName filesep 'optimization_parameters.csv'];
writematrix(Opt_param,storePath)
%/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/%