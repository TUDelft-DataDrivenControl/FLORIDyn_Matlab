%% Closed loop control main funtion
% Combines the EnKF with the MPC, both using FLORIDyn

pathToSimulation = ['2024_CLC' filesep '10T_HKN_06_to_09_221deg'];

%% Initialize
CLC_init;

%% Add UQ
%addpath(genpath("UQ_CLC"))

%% Load Reference data (temporary)

powSOWFA = importSOWFAFile([Sim.PathToSim 'Data' filesep ...
    'SOWFA_generatorPower.csv']);
powSOWFA(:,3) = powSOWFA(:,3)./(10^6*paramFLORIS.airDen);


uSOWFA = 0;
%% Run
% Run condition
run_controller  = true;

% "Global" time keeper
time.it         = 0;                % Current iteration
time.current    = Sim.StartTime;    % World time (LES + precursor time) 
time.current_nO = 0;                % World time reset to 0
time.end        = Sim.EndTime;      % World time to end controller
time.stepWorld  = 0.5;              % World / LES time step 
time.stepEnKF   = Sim.TimeStep * EnKF.nS;       % Time step of the EnKF
time.stepMPC    = Sim.TimeStep * CLC.Time.nS;   % Time step of the MPC

while run_controller
    % Increase counter & check break condition
    disp(['t = ' num2str(time.current,7)])
    time.it = time.it + 1;

    if time.current >= time.end; run_controller = false; end % Simulation end

    %% Estimator
    if mod(time.current_nO, time.stepEnKF) == 0 && time.current > Sim.StartTime
        disp('////////// Run EnKF //////////')
        % ======= Run Ensemble Kalman Filter ======= 
        CLC_EnKF_Par;
    end

    %% MPC
    if mod(time.current_nO, time.stepMPC) == 0 && time.current > Sim.StartTime
        % ======= Advance ensembles to current time step
        disp('Collect ensemble states to run controller')
        if mod(time.current_nO, time.stepEnKF) ~= 0
            error("Control not calculated following an EnKF step!")
        end


        % ======= Collect "true state" from ensembles
        [T,M] = EnKF_getMeanState(EnKF, T, paramFLORIDyn);
        
        if EnKF.StoreData
            folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
                filesep 'EnKF_States'];
            if ~exist(folder_name, "dir")
                mkdir(folder_name)
            end
            writematrix(T.States_WF(:,1)',[folder_name filesep...
                'vel'],'WriteMode','append')
            writematrix(T.States_WF(:,2)',[folder_name filesep...
                'dir'],'WriteMode','append')
            writematrix(T.States_Var_Vel',[folder_name filesep...
                'vel_var'],'WriteMode','append')
            writematrix(T.States_Var_Dir',[folder_name filesep...
                'dir_var'],'WriteMode','append')
        end
        

        % DEBUGGING START
        CLC.Time.StartTime = Sim.StartTime;
        % DEBUGGING END

        % ======= MPC ======= 
        disp('////////// Run Controller ////////// ')
        % ===== Predict Wind conditions
        Wind_Con = predict_wind(EnKF.Wind, time.current, CLC, Sim.TimeStep);
        [Con_Test, CLC_Test] = controller(T, Wind_Con, Sim, Con, Vis, paramFLORIDyn,...
            paramFLORIS,CLC, time.current);
        
        %  ======= Store data ======= 
        if or(EnKF.StoreData, CLC.Con.StoreData)
            folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
                filesep 'FlowField'];
            if ~exist(folder_name, "dir")
                mkdir(folder_name)
            end
            %   Generate flow field and mask
            fieldLims = Vis.FlowField.Lims;
            fieldRes  = Vis.FlowField.Res;
            xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
            yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);
            [X,Y] = meshgrid(xAx,yAx);

            % Get Mask
            W = getWeightsT(X(:), Y(:), T,...
                Sim.Dyn.Vel, Sim.TimeStep);
            W = sum(W,2);
            W = reshape(W,size(X));

            % Get flow field data
            U = GetFlowField(T,Wind,Vis,paramFLORIS);
            
            writematrix(U,[folder_name filesep...
                pad(num2str(time.current),7,"left",'0') '_U.csv'])
            writematrix(W,[folder_name filesep...
                pad(num2str(time.current),7,"left",'0') '_W.csv'])
            if ~exist([folder_name filesep 'X.csv'],"file")
                writematrix(X,[folder_name filesep 'X.csv'])
                writematrix(Y,[folder_name filesep 'Y.csv'])
            end
            writematrix([time.current T.States_OP(:,1)'],...
                [folder_name filesep 'OP_x.csv'],"WriteMode","append")
            writematrix([time.current T.States_OP(:,2)'],...
                [folder_name filesep 'OP_y.csv'],"WriteMode","append")

            clear fieldLims fieldRes xAx yAx X Y W U
        end
    end
    

    %% Time
    time.current    = time.current    + time.stepWorld;
    time.current_nO = time.current_nO + time.stepWorld;

    
end

disp('Simulation done!')

%% Store data
if EnKF.StoreData
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'EnKF_Measurements'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end
    writetable(cat(1,EnKF.M{:}),[folder_name filesep...
         'EnKF_measurements.csv'])
end

%% 
if and(EnKF.StoreData, CLC.Con.StoreData)
    % Generates plots from the outputs in .jpg and vector .pdf format
    postprocess_output([Sim.PathToSim 'Results' filesep CLC.CaseName])
end




%% Plot EnKF power
%EnKF_Vis = plotK_CombinedMeasurement(EnKF_Vis, EnKF, T, Sim, paramFLORIS, 5, 1:10);

%% Plot Flow field
%PlotFlowField(T,[],Wind,Sim,Vis,paramFLORIDyn,paramFLORIS);
%colormap(flipud(mako(100))); clim([2,10])
%% Calculate outputs
%e = error_EnKF(EnKF,T.nT,powSOWFA,uSOWFA,30600);