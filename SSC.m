%% Closed loop control main funtion
% Combines the EnKF with the MPC, both using FLORIDyn

pathToSimulation = ['10T_HKN_het_p10ramp_201deg'];

%% Initialize
CLC_init;

%% Add UQ
%addpath(genpath("UQ_CLC"))

%% Load Reference data (temporary)

powSOWFA = 0;

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


%% Setup zeroMQ server & init SOWFA
zmqServer = zeromqObj('/home/marcusbecker/OpenFOAM/marcusbecker-2.4.0/jeromq/jeromq-0.4.4-SNAPSHOT.jar',43660,3600,true);

% Load the yaw setpoint LUT and set-up a simple function
nTurbs = 10;

% Initial control settings
yawAngleArrayOut  	= 225.0*ones(1,nTurbs); % Initial settings
torqueArrayOut 		= 0.0 * ones(1,nTurbs); % Not used -- placeholder
pitchAngleArrayOut 	= 0.0 * ones(1,nTurbs);

dataSend = setupZmqSignal(torqueArrayOut,yawAngleArrayOut,pitchAngleArrayOut);
%dataSend = setupZmqSignal(yawAngleArrayOut,torqueArrayOut,pitchAngleArrayOut);

firstRun = true;
            
disp('Entering wind farm controller loop...');

while 1
    %% Receive information from SOWFA
    dataReceived = zmqServer.receive();
    currentTime  = dataReceived(1,1);
    measurementData = dataReceived(1,2:end);

    % Measurements: [genPower,rotSpeedF,azimuth,rotThrust,rotTorque,genTorque,nacYaw,bladePitch]
    %generatorPowerArray = measurementVector(1:8:end);
    generatorPowerArray = measurementData(1:8:end);
    rotorSpeedArray     = measurementData(2:8:end);
    azimuthAngleArray   = measurementData(3:8:end);
    rotorThrustArray    = measurementData(4:8:end);
    rotorTorqueArray    = measurementData(5:8:end);
    genTorqueArray      = measurementData(6:8:end);
    nacelleYawArray     = measurementData(7:8:end);
    bladePitchArray     = measurementData(8:8:end);

    
    
    % Update message string
    t = currentTime;
    time.current = currentTime;
    % Increase counter & check break condition
    disp(['t = ' num2str(time.current,7)])
    time.it = time.it + 1;

    %% Estimator
    if mod(time.current, time.stepEnKF) == 0 && ~firstRun  %time.current > Sim.StartTime
        disp('////////// Run EnKF //////////')
        % ======= Run Ensemble Kalman Filter ======= 
        CLC_EnKF_Par;
    end

    %% First entries arrive at ______.dt
    if firstRun; dt = rem(currentTime,10); firstRun = false; end
    
    %% MPC
    if mod(time.current, time.stepMPC) == 0 && time.current > CLC.Time.StartTime
        % ======= Advance ensembles to current time step
        disp('Collect ensemble states to run controller')
        if mod(time.current, time.stepEnKF) ~= 0
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
        
        % ======= MPC ======= 
        disp('////////// Run Controller ////////// ')
        % ===== Update current turbine orientation
        CLC.g0 = 270-nacelleYawArray';

        % ===== Predict Wind conditions
        Wind_Con    = predict_wind(EnKF.Wind, time.current, CLC, Sim.TimeStep);
        [Con, CLC]  = controller(T, Wind_Con, Sim, Con, Vis, paramFLORIDyn,...
            paramFLORIS,CLC, time.current); % Sim.StartTime);

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


    %% Set the control outputs
    % if and(t >= 20600, t<=20900)  
	% 	% Turbines stay orthogonal to uniform 60 deg wind dir change) 
	% 	yawAngleArrayOut = ones(1,nTurbs)*interp1([20600,20900],[255,195],t);
    % end
	yawAngleArrayOut = getYaw(Con.YawData,(1:T.nT)',currentTime);

    %% Prepare message for SOWFA
    disp([datestr(rem(now,1)) '__    Synthesizing message string.']);
    dataSend = setupZmqSignal(torqueArrayOut,yawAngleArrayOut,pitchAngleArrayOut);
    % dataSend = setupZmqSignal(yawAngleArrayOut,torqueArrayOut,pitchAngleArrayOut); %torque went to 0
    

    if currentTime >= (time.end - dt)
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
    end
    %% System
    % run true system every time step
    disp('Update system')
    % Send a message (control action) back to SOWFA
    zmqServer.send(dataSend);

    %% Time
    time.current_nO = time.current_nO + dt;
    
end

% Close connection
zmqServer.disconnect()

function [dataOut] = setupZmqSignal(torqueSignals,yawAngles,pitchAngles)
	dataOut = [];
    for i = 1:length(yawAngles)
        dataOut = [dataOut torqueSignals(i) yawAngles(i) pitchAngles(i)];
    end
end