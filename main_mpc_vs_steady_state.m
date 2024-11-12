%% MainFLORIDyn Center-Line model
% Improved FLORIDyn approach over the gaussian FLORIDyn model
% Link folder with simulation data
pathToSimulation = ['2024_CLC' filesep '2T_cost_function_study'];%'2021_9T_Data'; %'2023_NAWEA_2T_yaw_steering';%'2021_3T_Case_yaw_ConstWind';%['2024_CLC' filesep '10T_HKN_C1_201deg']; %

%% Load data from the simulation
% Reset the Matlab Path and load essential paths & the simulation path
addPaths;

% Get the settings for the wind field, visualization, controller and Sim.
[Wind, Vis, Sim, Con] = setup();

%Sim.EndTime = Overwrite_EndTime; % <- Reclict from NAWEA 2023 runs

% Add according functions to the search path
addFLORISPaths;
addFLORIDynPaths;

% Load linked data
turbProp        = turbineArrayProperties();
paramFLORIS     = parameterFLORIS();
paramFLORIDyn   = parameterFLORIDyn();

%% Preprocess loaded data
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);
clear turbProp

%% ====== Init simulation
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

%% ============ RUN SIMULATION ============
% T := Simulation state (OP states, Turbine states, wind field states(OPs))
% M := Measurements from the simulation (T_red, T_addedTI, T_Ueff, T_pow)
changes_dir = -.03:.01:.03;
yaw_sp      = -30:1:30;
% changes_dir = 0;
% yaw_sp = -2:2:2;

for iD = 1:length(changes_dir)
for iY = 1:length(yaw_sp)
    % Manipulate yaw and dir
    delta_t = abs(yaw_sp(iY))/.3;
    if delta_t == 0; delta_t = 10; end
    Con.YawData(1:2,2) = 270 - 16.2;
    Con.YawData(3,1:2) = [500 + delta_t, 270 - yaw_sp(iY)];
    Con.YawData(4,2)   = Con.YawData(3,2);
    Con.YawData(4,1)   = 9000;
    
    Wind.Dir = [0, 270;...
        500, 270;...
        1000,500*changes_dir(iD)+270;...
        9000,500*changes_dir(iD)+270];

    % Run sim
    [T_tmp,M,Vis,Mint] = ...
        FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

    farmpower = sum(reshape(M.("Power generated [MW]"),2,[]));
    writematrix([changes_dir(iD) yaw_sp(iY) farmpower(101:end)],...
        [Sim.PathToSim 'Results' filesep 'experiments_yawed.csv'],"WriteMode","append")
    disp(iD*iY/(length(changes_dir)*length(yaw_sp))*100)
end
end




%% ============ RUN SIMULATION ============
% T := Simulation state (OP states, Turbine states, wind field states(OPs))
% M := Measurements from the simulation (T_red, T_addedTI, T_Ueff, T_pow)
mn = .2; % deg / 5s

n_trajectories = 21;
n_steps        = 400;
t_trajectories = 0:5:n_steps*5;
trajectories = cumsum([ones(n_trajectories,1)*270, randn(n_trajectories, n_steps)*mn],2);

% figure
% hold on
% for it = 1:n_trajectories
%     plot(0:5:500,trajectories(it,:))
% end


yaw_sp      = -30:1:30;
% changes_dir = 0;
% yaw_sp = -2:2:2;

for iT = 1:n_trajectories
for iY = 1:length(yaw_sp)
    % Manipulate yaw and dir
    delta_t = abs(yaw_sp(iY))/.3;
    if delta_t == 0; delta_t = 10; end
    Con.YawData(1:2,2) = 270;
    Con.YawData(3,1:2) = [500 + delta_t, 270 - yaw_sp(iY)];
    Con.YawData(4,2)   = Con.YawData(3,2);
    Con.YawData(4,1)   = 9000;
    
    Wind.Dir = [0, 270;...
        500+t_trajectories',trajectories(iT,:)'];

    % Run sim
    [T_tmp,M,Vis,Mint] = ...
        FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

    farmpower = sum(reshape(M.("Power generated [MW]"),2,[]));
    writematrix([iT yaw_sp(iY) farmpower(101:end)],...
        [Sim.PathToSim 'Results' filesep 'experiments_random_walk.csv'],"WriteMode","append")
    
end
disp(iT/n_trajectories*100)
end




%% ============ RUN Update Steady State SIMULATION ============
% T := Simulation state (OP states, Turbine states, wind field states(OPs))
% M := Measurements from the simulation (T_red, T_addedTI, T_Ueff, T_pow)
changes_dir = -.02:.005:.02;
yaw_sp      = -30:2.5:30;
opt_start   = 1; % 1 or 0

deg2upd     = 2; % Threshold to update




for iD = 1:length(changes_dir)
    % deg / (deg/s) = s
    if changes_dir(iD) ~= 0
        t2upd = deg2upd / abs(changes_dir(iD));
    else
        t2upd = 500;
    end
    t2upd = round(t2upd / Sim.TimeStep);

    % Reset simulation
    ori_current = [270-15*opt_start, 270];
    dir_current = 270;
    %tim_current = 0;
    dir_old     = dir_current;
    delta_t_T1  = abs(dir_current - dir_old)/.3;
    if (dir_current - dir_old) == 0; delta_t_T1=1;end

    Con.YawData = ...
        [0, 270-15*opt_start,270;
        500, 270-15*opt_start,270;
        1000, 270-15*opt_start,270;
        9000, 270-15*opt_start,270];


    WindDir_changing = [0, dir_current;...
                500, dir_current;...
                1500,1000*changes_dir(iD)+dir_current;...
                9000,(9000-500)*changes_dir(iD)+dir_current];

    %windDir_steady   = 
    iCon_off      = 2; %(iS-1)*3 + 2;

% BEGIN "SIMULATION"
for iS = 1:ceil(100/t2upd)
start_time    = 500 + (iS-1)*t2upd*Sim.TimeStep;
farm_power    = zeros(length(yaw_sp),length(500:5:1500));


    % ////// Part 1 //////
    % Optimization assuming STEADY STATE CONDITIONS
for iY = 1:length(yaw_sp)
    % Manipulate yaw and dir
    delta_t = abs((dir_current - yaw_sp(iY)) - ori_current(1))/.3;
    if delta_t == 0; delta_t = 10; end
    
    % Build Con.YawData
    timesteps = unique([delta_t delta_t_T1]);
    if length(timesteps)==1
        Con.YawData(end,:)=[];
    end
    Con.YawData((1:length(timesteps)) + iCon_off,1) = timesteps + start_time;
    
    F = griddedInterpolant([0,delta_t],[ori_current(1) dir_old-yaw_sp(iY)],...
        "linear","nearest");
    Con.YawData((1:length(timesteps)) + iCon_off,2) = F(timesteps);
    
    F = griddedInterpolant([0,delta_t_T1],[ori_current(2), dir_current],...
        "linear","nearest");
    Con.YawData((1:length(timesteps)) + iCon_off,3) = F(timesteps);

    Con.YawData(length(timesteps) + 1 + iCon_off,:) = ...
        [9000, Con.YawData(2 + iCon_off,2:3)];


    % Steady wind dir
    Wind.Dir = [WindDir_changing(1:2,:);
        start_time dir_current;...
        9000,dir_current];
    if start_time == WindDir_changing(2,1);Wind.Dir(3,1) = start_time + 1; end
    
    % Run sim
    [T_tmp,M,Vis,Mint] = ...
        FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

    farmpower = sum(reshape(M.("Power generated [MW]"),2,[]));
    farm_power(iY,:) = farmpower(101:end);

    writematrix([iS changes_dir(iD) yaw_sp(iY) farmpower(101:end)],...
        [Sim.PathToSim 'Results' filesep 'experiments_yawed_' ...
        num2str(t2upd) '_steps_ddir_' num2str(changes_dir(iD)) '_steady_state.csv'],...
        "WriteMode","append")
end



% Determine optimal yaw according to steady state
test_yaw = linspace(yaw_sp(1), yaw_sp(end), 101);
%[~, I]   = max(interp1(yaw_sp, farm_power(:,end), test_yaw, "spline"));
[~, I]   = max(farm_power(:,end));
%opt_yaw  = test_yaw(I);
opt_yaw  = yaw_sp(I);
disp(['Optimal steady state yaw: ' num2str(opt_yaw) ' deg'])

% Write optimal yaw time series:
% [iS, delta t,  start orientation, end orientation]
delta_t_T0 = abs((dir_current - opt_yaw) - ori_current(1))/.3;
writematrix([ ...
    iS 0 delta_t_T0,       ori_current(1), dir_current - opt_yaw; ...
    iS 1 delta_t_T1,    ori_current(2)  dir_current], ...
    [Sim.PathToSim 'Results' filesep 'experiments_yawed_' ...
        num2str(t2upd) '_steps_ddir_' num2str(changes_dir(iD)) '_opt_yaw.csv'],...
        "WriteMode","append")



for iY = 1:length(yaw_sp)
    % Manipulate yaw and dir
    delta_t = abs((dir_current - yaw_sp(iY)) - ori_current(1))/.3;
    if delta_t == 0; delta_t = 10; end
    
    
    % Build Con.YawData
    timesteps = unique([delta_t delta_t_T1]);
    if length(timesteps)==1
        Con.YawData(end,:)=[];
    end
    Con.YawData((1:length(timesteps)) + iCon_off,1) = timesteps + start_time;
    
    F = griddedInterpolant([0,delta_t],[ori_current(1) dir_old-yaw_sp(iY)],...
        "linear","nearest");
    Con.YawData((1:length(timesteps)) + iCon_off,2) = F(timesteps);
    
    F = griddedInterpolant([0,delta_t_T1],[ori_current(2), dir_current],...
        "linear","nearest");
    Con.YawData((1:length(timesteps)) + iCon_off,3) = F(timesteps);

    Con.YawData(length(timesteps) + 1 + iCon_off,:) = ...
        [9000, Con.YawData(2 + iCon_off,2:3)];



    % Changing wind dir
    Wind.Dir = WindDir_changing;
    
    % Run sim
    [T_tmp,M,Vis,Mint] = ...
        FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

    farmpower = sum(reshape(M.("Power generated [MW]"),2,[]));
    writematrix([iS changes_dir(iD) yaw_sp(iY) farmpower((101:101+t2upd-1)+(iS-1)*t2upd)],...
        [Sim.PathToSim 'Results' filesep 'experiments_yawed_' ...
        num2str(t2upd) '_steps_ddir_' num2str(changes_dir(iD)) '.csv'],...
        "WriteMode","append")
end

% Save data for next step
dir_old     = interp1(WindDir_changing(:,1), WindDir_changing(:,2), ...
    (iS-1)*t2upd*Sim.TimeStep+500);
dir_current = interp1(WindDir_changing(:,1), WindDir_changing(:,2), ...
    iS*t2upd*Sim.TimeStep+500);
delta_t_T1  = abs(dir_current - dir_old)/.3;
if delta_t_T1 == 0; delta_t_T1=1;end

%   Con.YawData
timesteps = unique([delta_t_T0 delta_t_T1 delta_t_T1+t2upd*Sim.TimeStep t2upd*Sim.TimeStep]);
Con.YawData((1:length(timesteps)) + iCon_off,1) = timesteps + start_time;

F = griddedInterpolant([0,delta_t_T0],[ori_current(1) dir_old-opt_yaw],...
    "linear","nearest");
Con.YawData((1:length(timesteps)) + iCon_off,2) = F(timesteps);

F = griddedInterpolant([0,t2upd*Sim.TimeStep,delta_t_T1+t2upd*Sim.TimeStep],[dir_old, dir_old, dir_current],...
    "linear","nearest");
Con.YawData((1:length(timesteps)) + iCon_off,3) = F(timesteps);

% Update current orientation
ori_current = Con.YawData(length(timesteps)+iCon_off-1,2:3);  %[dir_old-opt_yaw, dir_old];

iCon_off = iCon_off + length(timesteps) - 1;

disp(['Finished step Nr ' num2str(iS) ' out of ' num2str(ceil(100/t2upd))])
end

writematrix(Con.YawData, [Sim.PathToSim 'Results' filesep 'experiments_yawed_' ...
        num2str(t2upd) '_steps_ddir_' num2str(changes_dir(iD)) '_ConYawData.csv'])

disp(['////// Finished iD=' num2str(iD) ' //////'])

end



disp('////////////// SWITCHING TO MPC //////////////////')





% ============ RUN Update MPC SIMULATION ============
% T := Simulation state (OP states, Turbine states, wind field states(OPs))
% M := Measurements from the simulation (T_red, T_addedTI, T_Ueff, T_pow)
changes_dir = -.02:.005:.02;
yaw_sp      = -30:2.5:30;


deg2upd     = 2; % Threshold to update

% Ph in s (needs to be divisble by 5 (time step))
ph = 350;


for iD = 1:length(changes_dir)
    % deg / (deg/s) = s
    if changes_dir(iD) ~= 0
        t2upd = deg2upd / abs(changes_dir(iD));
    else
        t2upd = 500;
    end
    t2upd = round(t2upd / Sim.TimeStep);

    % Reset simulation
    ori_current = [270-10*opt_start, 270];
    dir_current = 270;
    %tim_current = 0;
    dir_old     = dir_current;
    delta_t_T1  = abs(dir_current - dir_old)/.3;
    if (dir_current - dir_old) == 0; delta_t_T1=1;end

    Con.YawData = ...
        [0, 270-10*opt_start,270; 
        500, 270-10*opt_start,270; 
        1000, 270-10*opt_start,270; 
        9000, 270-10*opt_start,270];


    WindDir_changing = [0, dir_current;...
                500, dir_current;...
                1500,1000*changes_dir(iD)+dir_current;...
                9000,(9000-500)*changes_dir(iD)+dir_current];

    %windDir_steady   = 
    iCon_off      = 2; %(iS-1)*3 + 2;

% BEGIN "SIMULATION"
for iS = 1:ceil(100/t2upd)
start_time    = 500 + (iS-1)*t2upd*Sim.TimeStep;
farm_power    = zeros(length(yaw_sp),length(500:5:1500));


    % ////// Part 1 //////
    % Optimization assuming STEADY STATE CONDITIONS
for iY = 1:length(yaw_sp)
    % Manipulate yaw and dir
    delta_t = abs((dir_current - yaw_sp(iY)) - ori_current(1))/.3;
    if delta_t == 0; delta_t = 10; end
    

    % Build Con.YawData
    timesteps = unique([delta_t delta_t_T1]);
    if length(timesteps)==1
        Con.YawData(end,:)=[];
    end
    Con.YawData((1:length(timesteps)) + iCon_off,1) = timesteps + start_time;
    
    F = griddedInterpolant([0,delta_t],[ori_current(1) dir_old-yaw_sp(iY)],...
        "linear","nearest");
    Con.YawData((1:length(timesteps)) + iCon_off,2) = F(timesteps);
    
    F = griddedInterpolant([0,delta_t_T1],[ori_current(2), dir_current],...
        "linear","nearest");
    Con.YawData((1:length(timesteps)) + iCon_off,3) = F(timesteps);

    Con.YawData(length(timesteps) + 1 + iCon_off,:) = ...
        [9000, Con.YawData(2 + iCon_off,2:3)];


    % Steady wind dir
    Wind.Dir = [WindDir_changing(1:2,:);
        start_time dir_current;...
        9000,dir_current];
    if start_time == WindDir_changing(2,1);Wind.Dir(3,1) = start_time + 1; end
    
    % Run sim
    [T_tmp,M,Vis,Mint] = ...
        FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

    farmpower = sum(reshape(M.("Power generated [MW]"),2,[]));
    farm_power(iY,:) = farmpower(101:end);

    writematrix([iS changes_dir(iD) yaw_sp(iY) farmpower(101:end)],...
        [Sim.PathToSim 'Results' filesep 'experiments_yawed_' ...
        num2str(t2upd) '_steps_ddir_' num2str(changes_dir(iD)) '_ph=' num2str(ph) '_sssim.csv'],...
        "WriteMode","append")
end

farm_energy = cumsum(farm_power,2).*5;

% Determine optimal yaw according to steady state
%test_yaw = linspace(yaw_sp(1), yaw_sp(end), 101);
%[~, I]   = max(interp1(yaw_sp, farm_power(:,end), test_yaw, "spline"));
[~, I]   = max(farm_energy(:,(start_time-500)/5 + ph/5));
%opt_yaw  = test_yaw(I);
opt_yaw  = yaw_sp(I);
disp(['Optimal mpc yaw: ' num2str(opt_yaw) ' deg'])

% Write optimal yaw time series:
% [iS, delta t,  start orientation, end orientation]
delta_t_T0 = abs((dir_current - opt_yaw) - ori_current(1))/.3;
writematrix([ ...
    iS 0 delta_t_T0,       ori_current(1), dir_current - opt_yaw; ...
    iS 1 delta_t_T1,    ori_current(2)  dir_current], ...
    [Sim.PathToSim 'Results' filesep 'experiments_yawed_' ...
        num2str(t2upd) '_steps_ddir_' num2str(changes_dir(iD)) '_opt_yaw_ph=' num2str(ph) '.csv'],...
        "WriteMode","append")



for iY = 1:length(yaw_sp)
    % Manipulate yaw and dir
    delta_t = abs((dir_current - yaw_sp(iY)) - ori_current(1))/.3;
    if delta_t == 0; delta_t = 10; end
    

    % Build Con.YawData
    timesteps = unique([delta_t delta_t_T1]);
    if length(timesteps)==1
        Con.YawData(end,:)=[];
    end
    Con.YawData((1:length(timesteps)) + iCon_off,1) = timesteps + start_time;
    
    F = griddedInterpolant([0,delta_t],[ori_current(1) dir_old-yaw_sp(iY)],...
        "linear","nearest");
    Con.YawData((1:length(timesteps)) + iCon_off,2) = F(timesteps);
    
    F = griddedInterpolant([0,delta_t_T1],[ori_current(2), dir_current],...
        "linear","nearest");
    Con.YawData((1:length(timesteps)) + iCon_off,3) = F(timesteps);

    Con.YawData(length(timesteps) + 1 + iCon_off,:) = ...
        [9000, Con.YawData(2 + iCon_off,2:3)];


    % Changing wind dir
    Wind.Dir = WindDir_changing;
    
    % Run sim
    [T_tmp,M,Vis,Mint] = ...
        FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

    farmpower = sum(reshape(M.("Power generated [MW]"),2,[]));
    writematrix([iS changes_dir(iD) yaw_sp(iY) farmpower((101:101+t2upd-1)+(iS-1)*t2upd)],...
        [Sim.PathToSim 'Results' filesep 'experiments_yawed_' ...
        num2str(t2upd) '_steps_ddir_' num2str(changes_dir(iD)) '_ph=' num2str(ph) '_dynsim.csv'],...
        "WriteMode","append")
end

% Save data for next step
dir_old     = interp1(WindDir_changing(:,1), WindDir_changing(:,2), ...
    (iS-1)*t2upd*Sim.TimeStep+500);
dir_current = interp1(WindDir_changing(:,1), WindDir_changing(:,2), ...
    iS*t2upd*Sim.TimeStep+500);
delta_t_T1  = abs(dir_current - dir_old)/.3;
if delta_t_T1 == 0; delta_t_T1=1;end

%   Con.YawData
timesteps = unique([delta_t_T0 delta_t_T1 delta_t_T1+t2upd*Sim.TimeStep t2upd*Sim.TimeStep]);
Con.YawData((1:length(timesteps)) + iCon_off,1) = timesteps + start_time;

F = griddedInterpolant([0,delta_t_T0],[ori_current(1) dir_old-opt_yaw],...
    "linear","nearest");
Con.YawData((1:length(timesteps)) + iCon_off,2) = F(timesteps);

F = griddedInterpolant([0,t2upd*Sim.TimeStep,delta_t_T1+t2upd*Sim.TimeStep],[dir_old, dir_old, dir_current],...
    "linear","nearest");
Con.YawData((1:length(timesteps)) + iCon_off,3) = F(timesteps);

% Update current orientation
ori_current = Con.YawData(length(timesteps)+iCon_off-1,2:3);  %[dir_old-opt_yaw, dir_old];

iCon_off = iCon_off + length(timesteps) - 1;

disp(['Finished step Nr ' num2str(iS) ' out of ' num2str(ceil(100/t2upd))])
end

writematrix(Con.YawData, [Sim.PathToSim 'Results' filesep 'experiments_yawed_' ...
        num2str(t2upd) '_steps_ddir_' num2str(changes_dir(iD)) '_ph=' num2str(ph) '_ConYawData.csv'])

disp(['////// Finished iD=' num2str(iD) ' //////'])

end




%% Test against one another
cols      = RdBu(length(changes_dir));
cols(5,:) = [0,0,0];
nC        = size(cols,1);

changes_path = 'Simulations/2024_CLC/2T_cost_function_study/Results/';


changes_names_ss = {
    'experiments_yawed_20_steps_ddir_-0.02',
    'experiments_yawed_27_steps_ddir_-0.015',
    'experiments_yawed_40_steps_ddir_-0.01',
    'experiments_yawed_80_steps_ddir_-0.005',
    'experiments_yawed_100_steps_ddir_0'
    'experiments_yawed_80_steps_ddir_0.005',
    'experiments_yawed_40_steps_ddir_0.01',
    'experiments_yawed_27_steps_ddir_0.015',
    'experiments_yawed_20_steps_ddir_0.02'};

changes_names_ph360 = {
    'experiments_yawed_20_steps_ddir_-0.02_ph=350',
    'experiments_yawed_27_steps_ddir_-0.015_ph=350',
    'experiments_yawed_40_steps_ddir_-0.01_ph=350',
    'experiments_yawed_80_steps_ddir_-0.005_ph=350',
    'experiments_yawed_100_steps_ddir_0_ph=350'
    'experiments_yawed_80_steps_ddir_0.005_ph=350',
    'experiments_yawed_40_steps_ddir_0.01_ph=350',
    'experiments_yawed_27_steps_ddir_0.015_ph=350',
    'experiments_yawed_20_steps_ddir_0.02_ph=350'};

% Changing wind dir
Wind.Dir = WindDir_changing;

figure

hold on


for iD = 1:length(changes_dir)
Wind.Dir = [0, 270;...
                500, 270;...
                1500,1000*changes_dir(iD)+270;...
                9000,(9000-500)*changes_dir(iD)+270];

Con.YawData = ...
    readmatrix([changes_path changes_names_ss{iD} '_ConYawData.csv']);

% Run sim
[~,M_ss,Vis,~] = ...
    FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);


%Con.YawData = readmatrix("Simulations/2024_CLC/2T_cost_function_study/Results/experiments_yawed_20_steps_ddir_-0.02_ph=350_ConYawData.csv");
Con.YawData = readmatrix([changes_path changes_names_ph360{iD} '_ConYawData.csv']);
% Run sim
[~,M_mpc,Vis,~] = ...
    FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

Con.YawData(:,2) = Con.YawData(:,3);
% Run sim
[~,M_bl,Vis,~] = ...
    FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

% Pure energy
% figure
% plot(M_ss.("Time [s]")(1:2:end)-500, sum(reshape(M_ss.("Power generated [MW]"),2,[])), '-', 'LineWidth', 2, 'Color', [0.9686    0.7176    0.6000])
% hold on
% plot(M_mpc.("Time [s]")(1:2:end)-500, sum(reshape(M_mpc.("Power generated [MW]"),2,[])), '-', 'LineWidth', 2, 'Color', [0.6452    0.8095    0.8913])
% plot(M_bl.("Time [s]")(1:2:end)-500, sum(reshape(M_bl.("Power generated [MW]"),2,[])), ':k', 'LineWidth', 2)
% hold off
% grid on
% xlim([0,500])
% legend({'Steady State', 'MPC','BL'})
% ylabel('Farm Power in MW')
% xlabel('Time in s')

% Normalised

bl_energy = cumsum(sum(reshape(M_bl.("Power generated [MW]"),2,[])));
plot3(M_ss.("Time [s]")(1:2:end)-500, ones(size(bl_energy))*changes_dir(iD),...
    cumsum(sum(reshape(M_ss.("Power generated [MW]"),2,[])))./bl_energy, '-', ...
    'LineWidth', 1.5, 'Color', cols(iD,:))

plot3(M_mpc.("Time [s]")(1:2:end)-500, ones(size(bl_energy))*changes_dir(iD),...
    cumsum(sum(reshape(M_mpc.("Power generated [MW]"),2,[])))./bl_energy, ...
    '--', ...
    'LineWidth', 1.5, 'Color', cols(iD,:))

fill3([M_mpc.("Time [s]")(1:2:end)-500; flipud(M_mpc.("Time [s]")(1:2:end)-500)]', ...
    [ones(size(bl_energy))*changes_dir(iD), ones(size(bl_energy))*changes_dir(iD)], ...
    [cumsum(sum(reshape(M_ss.("Power generated [MW]"),2,[])))./bl_energy, ...
    fliplr(cumsum(sum(reshape(M_mpc.("Power generated [MW]"),2,[])))./bl_energy)], ...
    cols(iD,:),'FaceAlpha',.29,'EdgeColor','none')
% plot(M_bl.("Time [s]")(1:2:end)-500, cumsum(sum(reshape(M_bl.("Power generated [MW]"),2,[])))./bl_energy, ...
%     '-k', 'LineWidth', 1)


end
hold off

grid on
xlim([0,500])
legend({'Steady State', 'MPC','BL'})
zlabel('Farm Energy in MW/MW')
ylabel('Dir change in deg/s')
xlabel('Time in s')
title('- steady state, -- mpc (ph 350)')
%%
PlotMeasurements(T_tmp,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS);