%% FLORIDyn function for closed-loop-control with FLORIDyn
% Version 0.1
%   Features:
%       - Yaw control
%       - Uses FLORIDyn as optimization model
%       - Uses FLORIDyn as true model
%% ====== Choose simulation ======
%pathToSimulation = '2023_9T_OL_yaw_steering';
pathToSimulation = '2023_2T_yaw_steering_demo';
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

addCLCPaths;
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
CLC.x0 = rand(CLC.Con.DOFpT*T.nT,1);

Sim.nSimSteps = CLC.Time.nS;

Pwr_out = zeros(CLC.Time.Sections*CLC.Time.nS,T.nT+2);
Ori_out = zeros(CLC.Time.Sections*CLC.Time.nS,T.nT+1);

%%
Pwr_out(:,end) = sum(...
    reshape(M_idealBL.("Power generated [MW]")(1:end-T.nT),T.nT,[]),1)';

%% ====== Run simulation ======
for iCLC = 1:CLC.Time.Sections
    % ESTIMATION
    %   The from here on used state is considered the best estimate of the
    %   true system state
    
    % Optimize
    [Con, CLC] = controller(T, Wind, Sim, Con, Vis, paramFLORIDyn,...
        paramFLORIS,CLC, Sim.StartTime);

    % Run until the next control update
    [T, M, Vis] = FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
    
    % Update time
    Sim.StartTime = Sim.StartTime + CLC.Time.SecDur;

    % Store Power
    %   Timestamps
    Pwr_out((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,1) = ...
        M.("Time [s]")(1:T.nT:end);
    %   Power
    for iT = 1:T.nT
        Pwr_out((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,1+iT) = ...
            M.("Power generated [MW]")(iT:T.nT:end);
    end
    
    % Check if trajectories have been optimised, if not, skip plotting
    if Sim.StartTime<CLC.Time.StartTime+0.001; continue; end

    % Store yaw angle
    Ori_out((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,:) = ...
        Con.YawData(1:CLC.Time.nS,:);

    % Plot Results
    %% Power Generated
    figure(109)
    lay = getLayout(T.nT);
    for iT = 1:T.nT
        subplot(lay(1),lay(2),iT)
        hold on
        plot(M.("Time [s]")(iT:T.nT:end), ...
            M.("Power generated [MW]")(iT:T.nT:end),...
            'LineWidth',1.5,'Color',[108, 194, 74]./255)

        if Sim.StartTime < CLC.Time.StartTime + CLC.Time.SecDur + 0.1
            plot(M_idealBL.("Time [s]")(iT:T.nT:end), ...
                M_idealBL.("Power generated [MW]")(iT:T.nT:end),':k',...
                'LineWidth',1)
            xlabel('Time (s)')
            ylabel('Power Gen (MW)')
            title(['Turbine ' num2str(iT-1)])
            grid on
            xlim([CLC.Time.StartTime,CLC.Time.StartTime + ...
            CLC.Time.Sections * CLC.Time.SecDur])
        end
        hold off
    end

    %% Total power
    figure(110)
    hold on
    if Sim.StartTime < CLC.Time.StartTime + CLC.Time.SecDur + 0.1
        xlabel('Time (s)')
        ylabel('Total Power Gen (MW)')
        grid on
        xlim([CLC.Time.StartTime,CLC.Time.StartTime + ...
            CLC.Time.Sections * CLC.Time.SecDur])
    end

    timePow = Pwr_out((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,1);
    refPow = Pwr_out((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,end);
    genPow = sum(Pwr_out((1:CLC.Time.nS)+(iCLC-1)*CLC.Time.nS,2:end-1),2);
    plot(timePow(refPow<genPow),genPow(refPow<genPow),...
        '+','Color',[108, 194, 74]./255,'LineWidth',2)
    plot(timePow(refPow>genPow),genPow(refPow>genPow),...
        'o','Color',[224, 60, 49]./255,'LineWidth',2)
    plot(timePow,refPow,':k','LineWidth',1)
    hold off
    
    %% Trajectory
    figure(111)
    for iT = 1:T.nT
        subplot(lay(1),lay(2),iT)
        hold on
        if Sim.StartTime < CLC.Time.StartTime + CLC.Time.SecDur + 0.1
            xlabel('Time (s)')
            ylabel('Orientation (deg)')
            title(['Turbine ' num2str(iT-1)])
            grid on
            plot(Wind.Dir(:,1),Wind.Dir(:,2),':k');
            xlim([CLC.Time.StartTime,CLC.Time.StartTime + ...
            CLC.Time.Sections * CLC.Time.SecDur])
        end
        plot(Con.YawData(1:end-1,1),Con.YawData(1:end-1,iT+1),...
            'LineWidth',1.5,'Color',[0, 184, 200]./255)
        plot(Con.YawData(1:CLC.Time.nS+1,1),...
            Con.YawData(1:CLC.Time.nS+1,iT+1),...
            'LineWidth',2.5,'Color',[12, 35, 64]./255)
        hold off
    end
    pause(.1)
end
Pwr_out(:,end+1) = sum(Pwr_out(:,2:end-1),2);
