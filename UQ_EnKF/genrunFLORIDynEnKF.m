% This file is used to run and initialize all settings in FLORIDyn in order 
% to perform both sensitivity analysis and Bayesian calibration

% Output Y of FLORIDyn for input X
% X comprises of parameters defined in uq_sensitivity.m and uq_calibration
function Y = genrunFLORIDynEnKF(X)
%% MainFLORIDyn Center-Line model
% Improved FLORIDyn approach over the gaussian FLORIDyn model

%% Load data from the simulation
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

%% Assign X
% X1 velWeight_dw
Sim.Dyn.Vel.IterSigma_DW    = X(1);
% X2 velWeight_cw
Sim.Dyn.Vel.IterSigma_CW    = X(2);
% X3 velWeight_t
Sim.Dyn.Vel.IterSigma_time  = X(3);
% X4 dirWeight_dw
Sim.Dyn.Dir.IterSigma_DW    = X(4);
% X5 dirWeight_cw
Sim.Dyn.Dir.IterSigma_CW    = X(5);
% X6 dirWeight_t
Sim.Dyn.Dir.IterSigma_time  = X(6);
% X7 nE 
EnKF.nE = round(X(7));

path2Storage = 'O:\TUDelft\EnKF_Sensitivity_Data\';%'/Volumes/MarcusSSD/TUDelft/EnKF_Sensitivity_Data/';
pathToSimulation = '.\Simulations\2021_9T_Data\'; % './Simulations/2021_9T_Data/';
%%
global m sim_p
if mod(sim_p,5)==0
    Vis.CL.Store = true;
    Vis.CL.Run = num2str(sim_p);
else
    Vis.CL.Store = false;
end
%% Preprocess loaded data
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);
clear turbProp

%% ====== Init simulation
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

%% Process and init EnKF Data
[EnKF, Wind] = EnKF_PreProcessing(EnKF, T, Wind, Sim);

% Visualization Settings
Vis_Settings.FlowFieldError.Online  = Vis.FlowField.Error.Online;
Vis.FlowField.Error.Online          = 0;
EnKF_Vis = EnKF_Visualizer(Vis_Settings);


meanCDir = zeros(EnKF.Sim.Sections, T.nOP * T.nT * T.nT);
stddCDir = zeros(EnKF.Sim.Sections, T.nOP * T.nT * T.nT);
meanCVel = zeros(EnKF.Sim.Sections, T.nOP * T.nT * T.nT);
stddCVel = zeros(EnKF.Sim.Sections, T.nOP * T.nT * T.nT);

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

%% ============ RUN SIMULATION ============
% T := Simulation state (OP states, Turbine states, wind field states(OPs))
% M := Measurements from the simulation (T_red, T_addedTI, T_Ueff, T_pow)
%% ============ RUN SIMULATIONS ============
tic
for iS = 1:EnKF.Sim.Sections
    
    % Set inputs
    Sim.StartTime   = EnKF.Sim.StartTime + (iS-1) * EnKF.Sim.SecDur;
    Sim.EndTime     = EnKF.Sim.StartTime + iS * EnKF.Sim.SecDur;
    Sim.nSimSteps   = EnKF.nS;
    
    nE = EnKF.nE;
    EnKF_States_Vel = EnKF.States.Vel;
    EnKF_States_Dir = EnKF.States.Dir;
    EnKF_States_OP  = reshape(EnKF.States_OP,[],EnKF.nE);
    EnKF_M = EnKF.M;
    EnKF_Output_Pow = EnKF.Output.Pow;
    EnKF_Interaction = EnKF.Interaction;
    if isfield(T,'C_Vel')
        EnKF_weightedInteractionVel = EnKF.weightedInteractionVel;
    end
    if isfield(T,'C_Dir')
        EnKF_weightedInteractionDir = EnKF.weightedInteractionDir;
    end

    for iE = 1:nE
        % ===== Assign relevant ensemble states
        T_tmp = EnKF_AssignEnStates(EnKF, T, iE);
        
        % =========== Run simulation ================
        [T_tmp,M,~] = ...
            FLORIDynCL(T_tmp,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

        % ===== Stores relevant measurements and states
        %% Store updated wind field states
        EnKF_States_Vel(:,iE) = T_tmp.States_WF(:,1);
        EnKF_States_Dir(:,iE) = T_tmp.States_WF(:,2);

        %% Store OP states
        EnKF_States_OP(:,iE) = T_tmp.States_OP(:);

        %% Store measurements
        EnKF_M{iE} = [EnKF_M{iE}; M];
        EnKF_Output_Pow(:,iE) = ...
            table2array(EnKF_M{iE}(end-T.nT+1:end,6));

        %% Store Turbine interactions
%         for iT = 1:T_nT
%             EnKF_Interaction{iE,iT} = ...
%                 [T_tmp.dep{iT}',T_tmp.intOPs{iT},T_tmp.weight{iT}];
%         end
        if isfield(T_tmp,'C_Vel')
            EnKF_weightedInteractionVel{iE} = T_tmp.C_Vel;
        end
        if isfield(T_tmp,'C_Dir')
            EnKF_weightedInteractionDir{iE} = T_tmp.C_Dir;
        end
    end
    
    EnKF.States.Vel = EnKF_States_Vel;
    EnKF.States.Dir = EnKF_States_Dir;
    EnKF.States_OP  = reshape(EnKF_States_OP,[],EnKF.nStatesOP*EnKF.nE);
    EnKF.M = EnKF_M;
    EnKF.Output.Pow = EnKF_Output_Pow;
    EnKF.Interaction = EnKF_Interaction;
    EnKF.weightedInteractionVel = EnKF_weightedInteractionVel;
    EnKF.weightedInteractionDir = EnKF_weightedInteractionDir;

    %% Combination
    % Calculate "true" state & state error covariance
    if EnKF.Vel.Correct
        if projectOntoTrueState
            [EnKF,~] = EnKF_projectOntoTrueState(EnKF,Sim,T,1);
            [EnKF, truePos] = EnKF_projectOntoTrueState(EnKF,Sim,T,2);
            %EnKF.States_OP = repmat(truePos,1,EnKF.nE);
        end
        C_xx_Vel = EnKF_CalcStateCov(EnKF.States.Vel,EnKF.nE,T.StartI);
        
        [C_xy_Vel, C_yy_Vel] = ...
            EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE,...
             'Method', 'InputOutput', 'Output',EnKF.Output.Pow);
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
            % create localization covariance and multiply with state
            % error covariance matrix
            LocCov = GaspariAndCohn1999(EnKF.Vel.cutOffLength,distOPs);
            
            K = EnKF_CalcKalmanGain('cxy_cyy_cee_loc', ...
                'Cxy',C_xy_Vel,'Cyy',C_yy_Vel,'Cee',EnKF.Vel.C_ee_Vel,...
                'Loc',LocCov);

            EnKF.States.Vel(:,iE) = EnKF.States.Vel(:,iE) + ...
                K * (d_P - EnKF.Output.Pow(:,iE) + ...
                (randn(1,T.nT)*EnKF.Output.C_ee_Pow_Chol)');

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
                    'Loc',LocCov,'C',EnKF.weightedInteractionDir{iE});
                
                % Plot Localization influence
%                 if and(iS == 1, iE==1)
%                     EnKF_Vis = plotK_Localization(EnKF_Vis, T, Vis, ...
%                         EnKF.Dir.cutOffLength, 2);
%                 end
            else
                % No localization
                K = C_xx_Dir/(C*C_xx_Dir + EnKF.Dir.C_ee_Dir);
            end
 
            % Correction
%             EnKF_Vis = plotK_Correction_Quiver(EnKF_Vis, EnKF, T, Vis, K,...
%                 (d_Dir_j - EnKF.States.Dir(T.StartI,iE)), iE, iS, 2, Sim);
            EnKF.States.Dir(:,iE) = EnKF.States.Dir(:,iE) + ...
                K * (d_Dir_j - EnKF.States.Dir(T.StartI,iE));

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
    folderName = [path2Storage 'Runs' filesep ...
        pad(num2str(sim_p),5,'left','0') filesep 'C_Mat'];
    if not(isfolder(folderName)); mkdir(folderName); end

    data = cell2mat(EnKF.weightedInteractionVel);
    for iT = 1:9
        data_mean = mean(reshape(data(iT,:),[],EnKF.nE),2)';
        data_stdd = std(reshape(data(iT,:),[],EnKF.nE),[],2)';

        is = T.nOP * T.nT * (iT-1) + 1;
        ie = T.nOP * T.nT * iT;
        meanCDir(iS, is:ie) = data_mean;
        stddCDir(iS, is:ie) = data_stdd;

        writematrix(data_mean,...
            [folderName filesep 'C_Vel_mean_T' num2str(iT-1) '.csv'], ...
            "WriteMode","append")
        writematrix(data_stdd,...
            [folderName filesep 'C_Vel_stdd_T' num2str(iT-1) '.csv'], ...
            "WriteMode","append")
    end

    data = cell2mat(EnKF.weightedInteractionDir);
    for iT = 1:9
        data_mean = mean(reshape(data(iT,:),[],EnKF.nE),2)';
        data_stdd = std(reshape(data(iT,:),[],EnKF.nE),[],2)';
        
        is = T.nOP * T.nT * (iT-1) + 1;
        ie = T.nOP * T.nT * iT;
        meanCDir(iS, is:ie) = data_mean;
        stddCDir(iS, is:ie) = data_stdd;

        writematrix(data_mean,...
            [folderName filesep 'C_Dir_mean_T' num2str(iT-1) '.csv'], ...
            "WriteMode","append")
        writematrix(data_stdd,...
            [folderName filesep 'C_Dir_stdd_T' num2str(iT-1) '.csv'], ...
            "WriteMode","append")
    end
    
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
% %               EnKF_Vis = plotK_FlowField(EnKF_Vis, EnKF, T, Wind, Sim, ...
% %                   Vis, paramFLORIDyn, paramFLORIS, 1);
%             EnKF_Vis = plotK_WeightedFlowAndSOWFA(EnKF_Vis, EnKF, T,...
%                 Vis, Sim, paramFLORIS, Wind,d(1,2));
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
tmeas = toc;

%% Get Power
M_comb = cat(1,EnKF.M{:});
Data_raw = reshape(M_comb.("Power generated [MW]"),[],EnKF.nE);
Data_Mean = mean(Data_raw,2);
Data_Stdv = std(Data_raw,[],2);
%% Store X data
folderName = [path2Storage 'Runs' filesep ...
    pad(num2str(sim_p),5,'left','0') filesep 'Settings'];
if not(isfolder(folderName)); mkdir(folderName); end
writematrix(X,[folderName filesep 'X.csv'])

%% Store Power data
folderName = [path2Storage 'Runs' filesep ...
    pad(num2str(sim_p),5,'left','0') filesep 'Power'];
if not(isfolder(folderName)); mkdir(folderName); end
writematrix(Data_Mean,[folderName filesep 'PowerGen_mean.csv'])
writematrix(Data_Stdv,[folderName filesep 'PowerGen_stdv.csv'])
%% Set output
Y = [Data_Mean', Data_Stdv',...
    reshape(meanCDir,1,[]), reshape(stddCDir,1,[]), ...
    reshape(meanCVel,1,[]), reshape(stddCVel,1,[])];

sim_p = sim_p + 1;
message = [...
    'Progess: ' num2str(round(100*sim_p/m)) ' %'];%\n'...
    %'Remaining time: ' timeEst(tmeas,m-sim_p) ' (hh:mm:ss)'];
%dispstat(sprintf(message)) % Can be deleted
disp(message) % Can be deleted
end
