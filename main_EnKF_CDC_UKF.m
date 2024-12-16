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

clear all
%% Multi-Main for the FLORIDyn Center-Line model
% File to run the Enseble Kalman Filter as described in
%   [1] Data Assimilation - The ensemble kalman filter, Geir Evensen, 2nd
%   Edition, 2009, Springer, DOI: 10.1007/978-3-642-03711-5
% Link folder with simulation data
pathToSimulation = '2021_9T_Data';%'2022_9T_Data_EnKF_Reference';%
%% CDC code
%     ___  ____    ___
%    / __)(  _ \  / __)
%   ( (__  )(_) )( (__ 
%    \___)(____/  \___)
%% Load data from the simulation
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
UKF             = UKF_settings();

%% Preprocess loaded data
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);
clear turbProp

%% ====== Init simulation
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
%% Process and init EnKF Data
[UKF, Wind] = UKF_PreProcessing(UKF, T, Wind, Sim);

% Visualization Settings
Vis_Settings.FlowFieldError.Online  = Vis.FlowField.Error.Online;
Vis.FlowField.Error.Online          = 0;
EnKF_Vis = EnKF_Visualizer(Vis_Settings);

%% ============ LOADING IN VALIDATION DATA ==================
if isequal(pathToSimulation,'2022_9T_Data_EnKF_Reference')
    powSOWFA = readmatrix([Sim.PathToSim 'Data' filesep 'FLORIDyn_generatorPower.csv']);
else
    powSOWFA = importSOWFAFile([Sim.PathToSim 'Data' filesep ...
        'SOWFA_generatorPower.csv']);
    powSOWFA(:,3) = powSOWFA(:,3)./(10^6*paramFLORIS.airDen);
end

%% TRUE STATE PROJECTION
projectOntoTrueState = true;

%% ============ RUN SIMULATIONS ============
tStart = tic;
tSim = 0;
tCom = 0;
tCor = 0;
% Error stats
% 1: Average RMSE of the wind speed across the field
% 2: Average Error of the wind speed across the field
% 3: Std of the Error of the wind speed across the field
% 4: Percentage of wind speed errors within 1 std
% 5: Percentage of wind speed errors within 2 std
% 6: Percentage of wind speed errors within 3 std
% 7: Average OP position std
% 8: Average OP position std within wind farm bounds
% 9: Average OP position std outside wind farm bounds
err  = zeros(UKF.Sim.Sections,9);
for iS = 1:UKF.Sim.Sections
    tic
    % Set inputs
    Sim.StartTime   = UKF.Sim.StartTime + (iS-1) * UKF.Sim.SecDur;
    Sim.EndTime     = UKF.Sim.StartTime + iS * UKF.Sim.SecDur;
    Sim.nSimSteps   = UKF.nS;
    
    for iE = 1:UKF.nE
        % ===== Assign relevant ensemble states
        T = EnKF_AssignEnStates(UKF,T,iE);
        
        % =========== Run simulation ================
        [T,M,Vis] = ...
            FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
        
        % ===== Stores relevant measurements and states
        UKF = EnKF_StoreEnStates(UKF,M, T, iE);
    end
    tSim = tSim + toc;


    %% Combination
    tic
    % Calculate mean position and project states to mean position
    [UKF,~] = EnKF_projectOntoTrueState(UKF,Sim,T,1);
    [UKF, truePos] = EnKF_projectOntoTrueState(UKF,Sim,T,2);

    % Overwrite position
    UKF.States_OP = repmat(truePos,1,UKF.nE);
    
    % Distance of the OPs to one another
    distOPs = EnKF_distOPs(truePos(:,1:3),T.StartI);

    % Calc C_phi and C_u based on true states
    [C_u, C_phi] = EnKF_calcC_uC_phi(truePos,...
        mean(UKF.States.Dir,2), Sim.Dyn, T.posBase, T.nOP,...
        Sim.TimeStep);

    % Calculate Localization
    LocCov = GaspariAndCohn1999(UKF.Vel.cutOffLength,distOPs);

    % Get output
    %   Power based on projected states
    P = EnKF_calcPower(C_u, UKF, T.D, paramFLORIS);

    % Get measurements from system
    [d, UKF] = EnKF_GetWFInputs(UKF,Sim,T,paramFLORIS);
    d_P = interp1(powSOWFA(1:T.nT:end,2),reshape(powSOWFA(:,3),T.nT,[])',Sim.EndTime)';

    % State correlation
    %   Velocity & Power
    [C_xy_Vel, C_yy_Vel] = ...
        EnKF_CalcStateCov(UKF.States.Vel, UKF.nE,...
        'Method', 'InputOutput', 'Output',P);
    %   Wind direction
    [C_xy_Dir, C_yy_Dir] = ...
        EnKF_CalcStateCov(UKF.States.Dir, UKF.nE,...
        'Method', 'InputOutput', 'Output',UKF.States.Dir(T.StartI,:));

    K_pu = EnKF_CalcKalmanGain('cxy_cyy_cee_loc', ...
                'Cxy',C_xy_Vel,...          % State-to-Output cov. matrix
                'Cyy',C_yy_Vel,...          % Output covariance matrix
                'Cee',UKF.Vel.C_ee_Vel,... % Measurement covariance matrix
                'Loc',LocCov);              % Localization

    K_phi = EnKF_CalcKalmanGain('cxy_cyy_cee_loc', ...
                'Cxy',C_xy_Dir,...          % State-to-Output cov. matrix
                'Cyy',C_yy_Dir,...          % Output covariance matrix
                'Cee',UKF.Dir.C_ee_Dir,...  % Measurement covariance matrix
                'Loc',LocCov);              % Localization


    %% Correction & reseeding
    % Get std
    % THIS SHOULD INCLUDE WEIGHTING
    stdVel = std(UKF.States.Vel,[],2);
    stdDir = std(UKF.States.Dir,[],2);

    meaVel = mean(UKF.States.Vel,2);
    meaDir = mean(UKF.States.Dir,2);
    C_xx_vel = (meaVel - UKF.States.Vel) * (meaVel - UKF.States.Vel)'/(UKF.nE - 1);
    C_xx_dir = (meaDir - UKF.States.Dir) * (meaDir - UKF.States.Dir)'/(UKF.nE - 1);

    % New mean state
    UKF.States.Vel(:,iE) = mean(UKF.States.Vel,2) + ...
                K_pu * (d_P - mean(P,2));

    UKF.States.Dir(:,iE) = mean(UKF.States.Dir,2) + ...
                K_phi * (d(:,2) - mean(UKF.States.Dir(T.StartI,:),2));

    % Reseed
%     for iE = 1:UKF.L
%         if mod(iE,2)==0
%             UKF.States.Vel(:,iE+1) = UKF.States.Vel(:,iE+1) + stdVel;
%         else
%             UKF.States.Dir(:,iE+1) = UKF.States.Dir(:,iE+1) + stdDir;
%         end
%     end
% 
%     for iE = UKF.L+1:UKF.nE-1
%         if mod(iE,2)==0
%             UKF.States.Dir(:,iE+1) = UKF.States.Dir(:,iE+1) - stdDir;
%         else
%             UKF.States.Vel(:,iE+1) = UKF.States.Vel(:,iE+1) - stdVel;
%         end
%     end
    
    % Reseed randomly
%     for iE = 1:UKF.nE
%         randomError = randn(T.nOP*T.nT,1);
%         
%     end

    

    tCom = tCom + toc;
    
    %%
    disp(['Progress ' num2str(iS/UKF.Sim.Sections*100,4) ' %'])

    %% Store C data
%     folderName = [Sim.PathToSim 'Results' filesep 'C_Mat'];
%     if not(isfolder(folderName)); mkdir(folderName); end
%     writematrix(cell2mat(EnKF.weightedInteractionVel),...
%         [folderName filesep 'Vel.csv'], ...
%         'WriteMode','append')
%     writematrix(cell2mat(EnKF.weightedInteractionDir),...
%         [folderName filesep 'Dir.csv'], ...
%         'WriteMode','append')
%     folderName = [Sim.PathToSim 'Results' filesep 'C_Mat_T2_IFAC'];
%     if not(isfolder(folderName)); mkdir(folderName); end
%     writematrix(C_u(3,:),...
%         [folderName filesep 'Vel_T2.csv'], ...
%         'WriteMode','append')
%         writematrix(C_phi(3,:),...
%         [folderName filesep 'Dir_T2.csv'], ...
%         'WriteMode','append')
    %% Store states mean and stdd
%     folderName = [Sim.PathToSim 'Results' filesep 'states'];
%     if not(isfolder(folderName)); mkdir(folderName); end
% 
%     writematrix(mean(EnKF.States.Vel,2)',...
%         [folderName filesep 'Vel_mean.csv'], ...
%         'WriteMode','append')
%     writematrix(std(EnKF.States.Vel,[],2)',...
%         [folderName filesep 'Vel_stdd.csv'], ...
%         'WriteMode','append')
% 
%     writematrix(mean(EnKF.States.Dir,2)',...
%         [folderName filesep 'Dir_mean.csv'], ...
%         'WriteMode','append')
%     writematrix(std(EnKF.States.Dir,[],2)',...
%         [folderName filesep 'Dir_stdd.csv'], ...
%         'WriteMode','append')
% 
%     inBounds = and(truePos(:,1)<3000,truePos(:,2)<3000);
%     writematrix(inBounds',...
%         [folderName filesep 'in_bounds.csv'], ...
%         'WriteMode','append')
    %% Plot Flow Field
%     if iS> 44
%         EnKF_Vis = plotK_GaussianFlowField(EnKF_Vis, EnKF, T, Vis, Sim,1,iS);
%         
%     end
    %EnKF_Vis = plotK_WeightedFlowField(EnKF_Vis, EnKF, T, Vis, Sim,1,iS,8.2);
    %EnKF_Vis = plotK_WeightedFlowField(EnKF_Vis, EnKF, T, Vis, Sim,2,iS,d(1,2));

    %% Plot variance
%     EnKF_Vis = plotK_StateAndVariance(EnKF_Vis, EnKF, T, 2,d(1,2),Vis);
%     folderName = [Sim.PathToSim 'Results' filesep 'EnKF_Var' filesep 'Dir'];
%     if not(isfolder(folderName)); mkdir(folderName); end
%     saveas(gcf,[folderName filesep num2str(iS) '.png'])
% 
%     EnKF_Vis = plotK_StateAndVariance(EnKF_Vis, EnKF, T, 1, 8.2,Vis);
%     folderName = [Sim.PathToSim 'Results' filesep 'EnKF_Var' filesep 'Vel'];
%     if not(isfolder(folderName)); mkdir(folderName); end
%     saveas(gcf,[folderName filesep num2str(iS) '.png'])

%     EnKF_Vis = plotK_WeightedStateAndVariance(EnKF_Vis, EnKF, T, Vis,...
%         Sim, 1, 8.2);
%     saveas(gcf,[Sim.PathToSim 'Results' filesep 'EnKF_VarWeighted' filesep...
%         num2str(iS) '.png'])
    %% Forcing same starting positions for all ensembles
%     nS = length(T.Names_OP);
%     for iStates = 1:nS
%         EnKF.States_OP(:,iStates:nS:end) = ...
%             repmat(mean(EnKF.States_OP(:,iStates:nS:end),2),1,EnKF.nE);
%     end
    
    %% Plot Flow Field comparison to SOWFA
%     if EnKF_Vis.Settings.FlowFieldError.Online
%         if sum(Vis.FlowField.Error.Steps==Sim.EndTime)==1
% % %             EnKF_Vis = plotK_GaussianFlowFieldError(EnKF_Vis, EnKF, T,...
% % %                 Vis, Sim, iS, paramFLORIS, Wind);
%               [EnKF_Vis, err(iS,:)] = plotK_FlowField(EnKF_Vis, UKF, T, Wind, Sim, ...
%                   Vis, paramFLORIDyn, paramFLORIS, 1);
% %             EnKF_Vis = plotK_WeightedFlowAndSOWFA(EnKF_Vis, EnKF, T,...
% %                 Vis, Sim, paramFLORIS, Wind,d(1,2));
%         end
%     end


%     T.States_OP = truePos;
%     T.States_WF(:,1) = mean(EnKF.States.Vel,2);
%     T.States_WF(:,2) = mean(EnKF.States.Dir,2);
%     createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
%             Sim.EndTime,'Vel')
%     createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
%             Sim.EndTime,'Dir')
%     createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
%             Sim.EndTime,'Eff')
%     createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
%             Sim.EndTime,'OP')
end
t = toc(tStart);
timeErr = (1:length(err))*UKF.Sim.SecDur;
writematrix([timeErr', err],['Error_CDC_' num2str(UKF.nE) '.csv'])


disp(['Sec. per sim. step: ' num2str(t/(Sim.nSimSteps * ...
    UKF.Sim.Sections * UKF.nE)) ...
    ' with ' num2str(UKF.nE) ' ensembles, '...
    num2str(t/UKF.nE) ' s sim. time per ensemble,' ...
    ' total sim. time: ' num2str(t) ' s.'])
fprintf(['Simulation:  ' num2str(tSim) ' s\nCombination: ' num2str(tCom)...
    ' s\nCorrection:  ' num2str(tCor) ' s\n'])
EnKF_Vis = plotK_CombinedMeasurement(EnKF_Vis, UKF, T, Sim, paramFLORIS, 5, 1:9);
EnKF_Vis = plotK_CombinedMeasurement(EnKF_Vis, UKF, T, Sim, paramFLORIS, 1, 1:9);

% figure(9999)
% hold on
% plot(timeErr(err(:,1)>0),err(err(:,1)>0,1),'--*','LineWidth',2)
% grid on
% xlabel('Step (-)')
% ylabel('Averaged RMSE (ms^{-1})')
% hold off
% Notes
%   yaw currently forced to 0, issue is that the state is noisy,
%   application of the wind direction (and wind speed, but here less
%   relevant) should be also filtered by Lejeune filter. Also check if
%   angle in deg or rad, if in rad -> really bad, if in deg, fine.
%
%   Problem persists, not the same state evaluated. Idea: create surrogate
%   OPs at distinct downwind distances which represent the chaotic OPs and
%   are equal for all Turbines. Essentially a grid, more interpolation :(