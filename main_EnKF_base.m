%% Multi-Main for the FLORIDyn Center-Line model
% File to run the Enseble Kalman Filter as described in
%   [1] Data Assimilation - The ensemble kalman filter, Geir Evensen, 2nd
%   Edition, 2009, Springer, DOI: 10.1007/978-3-642-03711-5
% Link folder with simulation data
pathToSimulation = '2021_9T_Data';%'2022_9T_Data_EnKF_Reference';%
%% Base code (nothing deleted)
%    ____                
%   / __ )____ _________ 
%  / __  / __ `/ ___/ _ \
% / /_/ / /_/ (__  )  __/
%/_____/\__,_/____/\___/ 
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
EnKF            = EnKF_settings();

%% Preprocess loaded data
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);
clear turbProp

%% ====== Init simulation
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
%U = readmatrix('U.csv');
%% Process and init EnKF Data
[EnKF, Wind] = EnKF_PreProcessing(EnKF, T, Wind, Sim);

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
tic
for iS = 1:EnKF.Sim.Sections
    
    % Set inputs
    Sim.StartTime   = EnKF.Sim.StartTime + (iS-1) * EnKF.Sim.SecDur;
    Sim.EndTime     = EnKF.Sim.StartTime + iS * EnKF.Sim.SecDur;
    Sim.nSimSteps   = EnKF.nS;
    
    for iE = 1:EnKF.nE
        % ===== Assign relevant ensemble states
        T = EnKF_AssignEnStates(EnKF,T,iE);
        
        % =========== Run simulation ================
        [T,M,Vis] = ...
            FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
        
        % ===== Stores relevant measurements and states
        EnKF = EnKF_StoreEnStates(EnKF,M, T, iE);
    end
    


    %% Combination
    % Calculate "true" state & state error covariance
    if EnKF.Vel.Correct
        if projectOntoTrueState
            [EnKF,~] = EnKF_projectOntoTrueState(EnKF,Sim,T,1);
            [EnKF, truePos] = EnKF_projectOntoTrueState(EnKF,Sim,T,2);
            %EnKF.States_OP = repmat(truePos,1,EnKF.nE);

            % Calc C_phi and C_u based on true states
            [C_u, C_phi] = EnKF_calcC_uC_phi(truePos,...
                mean(EnKF.States.Dir,2), Sim.Dyn, T.posBase, T.nOP,...
                Sim.TimeStep);
        end
        C_xx_Vel = EnKF_CalcStateCov(EnKF.States.Vel,EnKF.nE,T.StartI);

        % Power based on projected states
        P = EnKF_calcPower(C_u, EnKF, T.D, paramFLORIS);
        
        %figure;plot(EnKF.Output.Pow(1,:));hold on;plot(P(1,:))
        [C_xy_Vel, C_yy_Vel] = ...
            EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE,...
             'Method', 'InputOutput', 'Output',P);%EnKF.Output.Pow);
    end
    if EnKF.Dir.Correct
        C_xx_Dir = EnKF_CalcStateCov(EnKF.States.Dir,EnKF.nE,T.StartI);
    end
    if EnKF.TI.Correct
        C_xx_TI = EnKF_CalcStateCov(EnKF.States.TI,EnKF.nE,T.StartI);
    end
    
    %% Get measurements from validation source
    %   Calls the functions to get wind speed, direction and amb.
    %   turbulence as the normal simulation would and gets [U,phi,TI]
    [d, EnKF] = EnKF_GetWFInputs(EnKF,Sim,T,paramFLORIS);

    r = EnKF_GetAveForeignReduction(EnKF.M,T.nT);
    d(:,1) = d(:,1)./(r*10^(-2));

    %d(:,1) = interp1(U(:,1),U(:,2:end),Sim.EndTime)';
    d_P = interp1(powSOWFA(1:T.nT:end,2),reshape(powSOWFA(:,3),T.nT,[])',Sim.EndTime)';

    %d(:,1) = d(:,1) + (randn(1,T.nT)*EnKF.Vel.C_ee_Vel_Chol)';
    %d(:,2) = d(:,2) + (randn(1,T.nT)*EnKF.Dir.C_ee_Dir_Chol)';
    % Add error
    %   [1] Eq. 4.35  d_j = d + e_j
    %   [1] Eq. 4.36  C_ee = ee^T;
    %   This should be a random walk model with the given covariance error
    
    %% Correct Ensemble states
    % [1] Eq.4.37 
    % x_a,j = x_f,j + C_f,xx * M^T (M C_f,xx M^T + C_ee)^-1 * (d_j-M*x_f,j)
    % x_a,j = x_f,j + K * (d_j-M*x_f,j)
    disp(['Progress ' num2str(iS/EnKF.Sim.Sections*100,4)])
    for iE = 1:EnKF.nE
        
        if sum([EnKF.Vel.loc, EnKF.Dir.loc, EnKF.TI.loc])>0
            OPs_tmp = ...
                EnKF.States_OP(:,EnKF.nStatesOP*(iE-1)+1:EnKF.nStatesOP*(iE-1)+3);
            distOPs = EnKF_distOPs(OPs_tmp,T.StartI);
        end
        
        if EnKF.Vel.Correct
            % Set measurements and pollute with C_ee_Vel
            d_Vel_j = d(:,1) + ...
                (randn(1,T.nT)*EnKF.Vel.C_ee_Vel_Chol)';
            
            
            % === ADAPTIVE C_xx_Vel CALCULATION WITH DIFFERENT C MATRIX ===
%             [C_xx_Vel, C] = EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE,...
%               'Method', 'AdaptiveC', 'Interaction', EnKF.Interaction, ...
%               'StartI', T.StartI, 'iE', iE, 'nT', T.nT, 'nOP', T.nOP);
            C_xx_Vel = EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE,...
              'Method', 'AdaptiveWeightedC',...
              'C', EnKF.weightedInteractionVel{iE});
            % =========================================================== %
            
            if EnKF.Vel.loc
                % create localization covariance and multiply with state
                % error covariance matrix
                LocCov = GaspariAndCohn1999(EnKF.Vel.cutOffLength,distOPs);
                
%                 K = EnKF_CalcKalmanGain('starti_cxx_cee_loc', ...
%                     'Cxx',C_xx_Vel,'Cee',EnKF.Vel.C_ee_Vel,...
%                     'Loc',LocCov,'StartI',T.StartI);
                
%                 K = EnKF_CalcKalmanGain('c_cxx_cee_loc', ...
%                     'Cxx',C_xx_Vel,'Cee',EnKF.Vel.C_ee_Vel,...
%                     'Loc',LocCov,'C',C);
                
                K = EnKF_CalcKalmanGain('c_cxx_cee_loc', ...
                    'Cxx',C_xx_Vel,'Cee',EnKF.Vel.C_ee_Vel,...
                    'Loc',LocCov,'C',EnKF.weightedInteractionVel{iE});
%                 % Plot Localization influence
%                 if and(iS == 1, iE==1)
%                     EnKF_Vis = plotK_Localization(EnKF_Vis, T, Vis, ...
%                         EnKF.Vel.cutOffLength, 1);
%                 end
            else
                % No localization
                K = EnKF_CalcKalmanGain('c_cxx_cee', ...
                    'Cxx',C_xx_Vel,'Cee',EnKF.Vel.C_ee_Vel,...
                    'C',EnKF.weightedInteractionVel{iE});
            end
            
            % Plotting of the combined correction
%             EnKF_Vis = plotK_combined(EnKF_Vis, T.nT, EnKF, Vis, K,...
%                 (d_Vel_j - EnKF.States.Vel(T.StartI,iE)),iE,1);
            
            % Correction
%             EnKF.States.Vel(:,iE) = EnKF.States.Vel(:,iE) + ...
%                 K * (d_Vel_j - EnKF.States.Vel(T.StartI,iE));
            
            K = EnKF_CalcKalmanGain('cxy_cyy_cee_loc', ...
                'Cxy',C_xy_Vel,'Cyy',C_yy_Vel,'Cee',EnKF.Vel.C_ee_Vel,...
                'Loc',LocCov);
%             EnKF_Vis = ...
%                 plotK_Correction(EnKF_Vis, EnKF, T, Vis, K,...
%                 (d_Vel_j - EnKF.States.Vel(T.StartI,iE)), iE, iS, 1);
            
            EnKF.States.Vel(:,iE) = EnKF.States.Vel(:,iE) + ...
                K * (d_P - EnKF.Output.Pow(:,iE) + ...
                (randn(1,T.nT)*EnKF.Output.C_ee_Pow_Chol)');

%             if projectOntoTrueState
%                 EnKF = EnKF_projectOntoEnsembleState(EnKF,1);
%             end
        end
        if EnKF.Dir.Correct
            % Set measurements
            d_Dir_j = d(:,2) + ...
                (randn(1,T.nT)*EnKF.Dir.C_ee_Dir_Chol)';
            
            if EnKF.Dir.loc
                % create localization covariance and multiply with state
                % error covariance matrix
                LocCov = GaspariAndCohn1999(EnKF.Dir.cutOffLength,distOPs);
                K = EnKF_CalcKalmanGain('c_cxx_cee_loc', ...
                    'Cxx',C_xx_Dir,'Cee',EnKF.Dir.C_ee_Dir,...
                    'Loc',LocCov,'C',C_phi);%EnKF.weightedInteractionDir{iE});
                
                % Plot Localization influence
%                 if and(iS == 1, iE==1)
%                     EnKF_Vis = plotK_Localization(EnKF_Vis, T, Vis, ...
%                         EnKF.Dir.cutOffLength, 2);
%                 end
            else
                % No localization
                K = C_xx_Dir/(C*C_xx_Dir + EnKF.Dir.C_ee_Dir);
            end
            
            % Plotting of the combined correction
%             EnKF_Vis = plotK_combined(EnKF_Vis, T.nT, EnKF, Vis, K,...
%                 (d_Dir_j - EnKF.States.Dir(T.StartI,iE)),iE,2);
%             EnKF_Vis = plotK_Correction(EnKF_Vis, EnKF, T, Vis, K,...
%                 (d_Dir_j - EnKF.States.Dir(T.StartI,iE)), iE, iS, 2);
            % Correction
%             EnKF_Vis = plotK_Correction_Quiver(EnKF_Vis, EnKF, T, Vis, K,...
%                 (d_Dir_j - EnKF.States.Dir(T.StartI,iE)), iE, iS, 2, Sim);
            EnKF.States.Dir(:,iE) = EnKF.States.Dir(:,iE) + ...
                K * (d_Dir_j - EnKF.States.Dir(T.StartI,iE));

%             if projectOntoTrueState
%                 EnKF = EnKF_projectOntoEnsembleState(EnKF,2);
%             end
        end
        if EnKF.TI.Correct
            % Set measurements
            d_TI_j = d(:,3) + ...
                (randn(1,T.nT)*EnKF.TI.C_ee_TI_Chol)';
            
            if EnKF.TI.loc
                % create localization covariance and multiply with state
                % error covariance matrix
                LocCov = GaspariAndCohn1999(EnKF.TI.cutOffLength,distOPs);
                
                K = (LocCov.*C_xx_TI)/(LocCov(T.StartI,:).*...
                    C_xx_TI(T.StartI,:) + EnKF.TI.C_ee_TI);
                
                % Plot Localization influence
                if and(iS == 1, iE==1)
                    EnKF_Vis = plotK_Localization(EnKF_Vis, T, Vis, ...
                        EnKF.TI.cutOffLength, 2);
                end
            else
                % No localization
                K = C_xx_TI/(C_xx_TI(T.StartI,:) + EnKF.TI.C_ee_TI);
            end
            
            % Plotting of the combined correction
            EnKF_Vis = plotK_combined(EnKF_Vis, T.nT, EnKF, Vis, K,...
                (d_TI_j - EnKF.States.TI(T.StartI,iE)),iE,3);
            
            % Correction
            EnKF.States.TI(:,iE) = EnKF.States.TI(:,iE) + ...
                K * (d_TI_j - EnKF.States.TI(T.StartI,iE));
        end
        
    end
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
    if EnKF_Vis.Settings.FlowFieldError.Online
        if sum(Vis.FlowField.Error.Steps==Sim.EndTime)==1
% %             EnKF_Vis = plotK_GaussianFlowFieldError(EnKF_Vis, EnKF, T,...
% %                 Vis, Sim, iS, paramFLORIS, Wind);
              EnKF_Vis = plotK_FlowField(EnKF_Vis, EnKF, T, Wind, Sim, ...
                  Vis, paramFLORIDyn, paramFLORIS, 1);
%             EnKF_Vis = plotK_WeightedFlowAndSOWFA(EnKF_Vis, EnKF, T,...
%                 Vis, Sim, paramFLORIS, Wind,d(1,2));
        end
    end


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
t = toc;
disp(['Sec. per sim. step: ' num2str(t/(Sim.nSimSteps * ...
    EnKF.Sim.Sections * EnKF.nE)) ...
    ' with ' num2str(EnKF.nE) ' ensembles, '...
    num2str(t/EnKF.nE) ' s sim. time per ensemble,' ...
    ' total sim. time: ' num2str(t) ' s.'])

EnKF_Vis = plotK_CombinedMeasurement(EnKF_Vis, EnKF, T, Sim, paramFLORIS, 5, 1:9);
EnKF_Vis = plotK_CombinedMeasurement(EnKF_Vis, EnKF, T, Sim, paramFLORIS, 1, 1:9);

% Notes
%   yaw currently forced to 0, issue is that the state is noisy,
%   application of the wind direction (and wind speed, but here less
%   relevant) should be also filtered by Lejeune filter. Also check if
%   angle in deg or rad, if in rad -> really bad, if in deg, fine.
%
%   Problem persists, not the same state evaluated. Idea: create surrogate
%   OPs at distinct downwind distances which represent the chaotic OPs and
%   are equal for all Turbines. Essentially a grid, more interpolation :(