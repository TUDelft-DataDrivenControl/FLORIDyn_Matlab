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

% This file is used to run and initialize all settings in FLORIDyn in order 
% to perform both sensitivity analysis and Bayesian calibration

% Output Y of FLORIDyn for input X
% X comprises of parameters defined in uq_sensitivity.m and uq_calibration
function Y = genrunFLORIDynCLC(X)
%% Closed loop control

%% CLC init
% Get the settings for the wind field, visualization, controller and Sim.
[Wind, Vis, Sim, Con] = setup();

% ====== Add according functions to the search path
addFLORISPaths;
addFLORIDynPaths;
addEnKFPaths;

% ====== Load linked data
turbProp        = turbineArrayProperties();
paramFLORIS     = parameterFLORIS();
paramFLORIDyn   = parameterFLORIDyn();
CLC             = clc_settings(Sim);
EnKF            = EnKF_settings();
addCLCPaths;

%% Assign X
% ========= FLORIDyn =========
% X(1)  = 'sig_vel_dw'; % D
Sim.Dyn.Vel.IterSigma_DW    = X(1) * 178.4;
disp(['In: Vel.IterSigma_DW = ' num2str(X(1))])
% X(2)  = 'sig_vel_cw'; % D
Sim.Dyn.Vel.IterSigma_CW    = X(2) * 178.4;
disp(['In: Vel.IterSigma_CW = ' num2str(X(2))])
% X(3)  = 'sig_vel_t';  % s
Sim.Dyn.Vel.IterSigma_time  = X(3);
disp(['In: Vel.IterSigma_time = ' num2str(X(3))])
% X(4)  = 'd'; 
Sim.Dyn.Advection           = X(4);
disp(['In: Advection = ' num2str(X(4))])

% ========= FLORIS ===========
% X(5)  = 'eta';    
% paramFLORIS.eta     = X(5);              % Tuned by hand 
% disp(['In: Eta = '  num2str(X(5))])
% X(6)  = 'pp'; 
paramFLORIS.p_p     = X(5);
disp(['In: Pp = '   num2str(X(5))])
% X(7)  = 'alpha'; 
paramFLORIS.alpha   = X(6);
disp(['In: alpha = ' num2str(X(6))])
% X(8)  = 'beta';
paramFLORIS.beta    = X(7);
disp(['In: beta = ' num2str(X(7))])
% X(9)  = 'ka';
paramFLORIS.k_a     = X(8);
disp(['In: k_a = '  num2str(X(8))])
% X(10) = 'kb'; 
% paramFLORIS.k_b     = X(10);              % Removed
% disp(['In: k_b = '  num2str(X(10))])
% % X(11) = 'kfa';                          % Removed
% paramFLORIS.k_fa    = X(11);
% disp(['In: k_fa = ' num2str(X(11))])
% X(12) = 'kfb';
paramFLORIS.k_fb    = X(9);
disp(['In: k_fb = ' num2str(X(9))])
% X(13) = 'kfc'; 
paramFLORIS.k_fc    = X(10);
disp(['In: k_fc = ' num2str(X(10))])
% X(14) = 'kfd'; 
% paramFLORIS.k_fd    = X(14);              % Removed
% disp(['In: k_fd = ' num2str(X(14))])
% X(15) = 'kTI'; 
paramFLORIS.TIexp   = X(11);
disp(['In: kTI = '  num2str(X(11))])
% ========= EnKF ===========
% X(16) = 'loc_vel'; % D     
EnKF.Vel.cutOffLength = X(12) * 178.4;
disp(['In: Vel.cutOffLength = ' num2str(X(12))])
% X(17) = 'pn_vel';  % m/s     
EnKF.Vel.ProcessNoise = X(13);
disp(['In: Vel.ProcessNoise = ' num2str(X(13))])
% X(18) = 'mn_pow';  % MW
EnKF.Output.PowNoiseVar = X(14);
disp(['In: PowNoiseVar = '      num2str(X(14))])

paramFLORIS.airDen  = 1.225;
%% Store inputs
if exist('UQ_Inputs.csv','file')
    % Add new inputs to inputs table
    writematrix(X,'UQ_Inputs.csv','WriteMode','append')
else
    % create files
    writematrix(X,'UQ_Inputs.csv')
end

%% 
% ====== Preprocess loaded data
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);

% ====== Init simulation
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

% ====== Process and init EnKF Data
[EnKF, Wind] = EnKF_PreProcessing(EnKF, T, Wind, Sim);

% ====== Visualization Settings
Vis_Settings.FlowFieldError.Online  = Vis.FlowField.Error.Online;
Vis.FlowField.Error.Online          = 0;
EnKF_Vis = EnKF_Visualizer(Vis_Settings);

% ====== Init Control parameters ======
CLC.gd = zeros(T.nT,1);
CLC.g0 = T.States_WF(T.StartI,2); % TODO Yaw_0 is WF & yaw
CLC.x0 = rand(CLC.Con.DOFpT*T.nT,1);

clear Vis_Settings turbProp


%% Global variables
% sim_p : Current simulation iteration
% m     : Overall number of simulations
global m sim_p


%% ============ RUN SIMULATION ============

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
    %disp(['t = ' num2str(time.current,7)])
    time.it = time.it + 1;

    if time.current >= time.end; run_controller = false; end % Simulation end

    %% Estimator
    if mod(time.current_nO, time.stepEnKF) == 0 && time.current > Sim.StartTime
        %disp('Run EnKF')
        % ======= Run Ensemble Kalman Filter ======= 
        %disp(['EnKF running, time: ' num2str(time.current) ' s'])
        CLC_EnKF_Par;

        % ======= Store data ======= 

    end

    %% MPC
    if mod(time.current_nO, time.stepMPC) == 0 && time.current > Sim.StartTime
        % ======= Advance ensembles to current time step
        %disp('Collect ensemble states to run controller')
        if mod(time.current_nO, time.stepEnKF) ~= 0
            error("Control not calculated following an EnKF step!")
        end
        % ======= Collect "true state" from ensembles
        T = EnKF_getMeanState(EnKF, T, paramFLORIDyn);
        
        % ======= MPC ======= 
        %disp('Run Controller')
        % [Con_Test, CLC_Test] = controller(T, Wind, Sim, Con, Vis, paramFLORIDyn,...
        %     paramFLORIS,CLC, Sim.StartTime);

        % ======= Store data ======= 
        
        % %   Store flow field with mask
        % fieldLims = Vis.FlowField.Lims;
        % fieldRes  = Vis.FlowField.Res;
        % xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
        % yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);
        % [X,Y] = meshgrid(xAx,yAx);
        % W = getWeightsT(X(:), Y(:), T,...
        %     Sim.Dyn.Vel, Sim.TimeStep);
        % W = sum(W,2);
        % W = reshape(W,size(X));
        % PlotFlowField(T,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS);
        % hold on
        % EnKF_plot_mask(X,Y,-W,-2,[1,1,1])
        % hold off
        % exportgraphics(gcf,[Sim.PathToSim filesep 'Results' ...
        %     filesep 'EnKF201deg' filesep 'FF_' num2str(time.current) '.png'])
        % close(gcf)
    end
    

    %% Time
    time.current    = time.current    + time.stepWorld;
    time.current_nO = time.current_nO + time.stepWorld;

    %% System
    % run true system every time step
    %disp('Update system')
    
end

sim_p = sim_p + 1;
message = [...
    'Progess: ' num2str(round(100*sim_p/m)) ' %%\n'];%...
    %'Remaining time: ' timeEst(tmeas,m-sim_p) ' (hh:mm:ss)'];
%dispstat(sprintf(message)) % Can be deleted
disp(message)

%% Get QoI
Y = error_EnKF(EnKF,T.nT,powSOWFA,uSOWFA,30600);


%% Store outputs
if exist('UQ_Outputs.csv','file')
    % Add new inputs to inputs table
    writematrix(Y','UQ_Outputs.csv','WriteMode','append')
else
    % create files
    writematrix(Y','UQ_Outputs.csv')
end
end
