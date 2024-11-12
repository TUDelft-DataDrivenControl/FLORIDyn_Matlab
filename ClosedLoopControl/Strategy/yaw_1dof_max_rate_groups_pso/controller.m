function [Con, CLC] = controller(T, Wind, Sim, Con, Vis, paramFLORIDyn, paramFLORIS,CLC, SimTime)
%CONTROLLER determines the control inputs for the near future and beyond.
%   It utilizes a given cost-function to evaluate the 
%% Bezier_yaw
% This controller is based on Bezier-Curves, which allow a smooth
% transition between distinct setpoints. These setpoints are subject to the
% optimization.
% Bezier-curves further allow the user to calculate the derivative over
% time, which allows us to enforce a yaw-rate limit (or axial induction)

%% Check if controller should be active
if SimTime<CLC.Time.StartTime
    return 
end

%% Task 1 
% Calculate the normalised rate limit
%   The rate limit for the normed space is calculated by multiplying the
%   rad/s with the time window, which leads to rad/[-].
rate_lim = CLC.Set.yawRateLimit * Sim.TimeStep * CLC.Con.horizon_action;
%   In future versions of the code, one could also change the action time
%   for each turbine, allowing for different action horizons.

%% Generate input for turbines - Step 1 
% Read starting point based on random input or previous step
x0 = CLC.x0;

% Control Update ratio
cu = CLC.Time.nS/CLC.Con.horizon_action;

% Constants
% Time line (+1 to account for current time step)
t_n = linspace(0,1,CLC.Con.horizon_action+1);
t   = linspace(SimTime, SimTime + CLC.Con.horizon_action * Sim.TimeStep,...
    CLC.Con.horizon_action+1);


[~,~,~,~,~,~,t2g,~] =...
    calcT(T, mean(T.States_WF(T.StartI,1)), CLC.Con.horizon_action, Sim.TimeStep);

t2g_time = zeros(size(t2g(:,1)));


for iG = 1:size(t2g,1)
    tic
    % Go through the different turbine groups that affect each other

    % Find out which optimization parameters are relevant to this group
    ind = reshape((1:T.nT*CLC.Con.DOFpT)',CLC.Con.DOFpT,[]);
    ind = ind(:,t2g(iG,:)>0);
    
    % Reduce wind farm to relevant turbines and data
    [T_tmp, Con_tmp, CLC_tmp] = reduceT(T, Con, CLC, t2g(iG,:)>0);
    CLC_tmp.x0 = CLC.x0(ind(:));
    CLC_tmp.g0 = CLC.g0(ind(:));
    
    disp(['Optimizing ' num2str(T_tmp.nT) ' turbine farm.'])
    % Build cost function wrapper that incorporates the reduced x into x
    f = @(x) cost_function_wrapper(x,T_tmp,Wind,Sim,Con_tmp,Vis,paramFLORIDyn,...
    paramFLORIS,CLC_tmp,rate_lim,SimTime,t_n,t,cu);

    % ================== Optimise the trajectories ========================== %
    options = optimoptions('particleswarm',...
        'SwarmSize',100,...
        'MaxIterations',CLC.max_iterations,...
        'MaxStallIterations',4,...
        'UseParallel',true,...
        'Display','iter');
    
    x = particleswarm(f,length(CLC_tmp.x0),...
            zeros(length(CLC_tmp.x0),1),...
            ones(length(CLC_tmp.x0),1),options);
    %x = rand(size(CLC_tmp.x0));

    CLC.x0(ind(:)) = x';

    t2g_time(iG) = toc;
    disp(['The optimization took ' num2str(t2g_time(iG)) 's.'])
end


opt_time = sum(t2g_time);
disp(['The optimization took ' num2str(opt_time) 's total.'])

%% Store optimised trajectory in Con.YawData
% Step 1 Generate the rajectory
[tr, CLC.g0] = mr_1_tr(t_n, CLC.x0, cu, rate_lim, CLC.g0);
Con.YawData = [t', tr;
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];


%% Store data
if CLC.Con.StoreData
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'Con_Trajectories'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end
    % Trajectories
    writematrix(Con.YawData,...
        [folder_name filesep pad(num2str(SimTime),7,"left",'0') '.csv'])
    
    % Convexity / Landscape
    % folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
    %     filesep 'Con_Landscape'];
    % if ~exist(folder_name, "dir")
    %     mkdir(folder_name)
    % end
    % f = @(x) cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    %     paramFLORIS,CLC,rate_lim,SimTime,t_n,t,cu);
    % 
    % J_values = landscape(CLC.x0,f);
    % writematrix(squeeze(J_values(:,:,1)),...
    %     [folder_name filesep pad(num2str(SimTime),7,"left",'0') '_J.csv'])
    % writematrix(squeeze(J_values(:,:,2)),...
    %     [folder_name filesep pad(num2str(SimTime),7,"left",'0') '_E1.csv'])
    % writematrix(squeeze(J_values(:,:,3)),...
    %     [folder_name filesep pad(num2str(SimTime),7,"left",'0') '_E2.csv'])
    
    % Control times
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'Con_Optimization_time'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end
    writematrix([SimTime, opt_time],...
        [folder_name filesep 'Opti_Time.csv'],"WriteMode","append")
    
    % Combined groups
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'Con_Optimization_groups'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end
    writematrix([t2g,t2g_time],...
        [folder_name filesep pad(num2str(SimTime),7,"left",'0') ...
        '_OptiGroups.csv'])

    % Measurements from optimal case
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'Con_Optim_Predictions'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end
    Sim_tmp = Sim;
    Sim_tmp.StartTime = SimTime;
    Sim_tmp.nSimSteps = CLC.Con.horizon_prediction;

    [~,M] = cost_function(T,Wind,Sim_tmp,Con,Vis,paramFLORIDyn,paramFLORIS);
    writetable(M,...
        [folder_name filesep pad(num2str(SimTime),7,"left",'0') '_M.csv'])
end
end











function J = cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rn,SimTime,t_n,t,cu)

% Correct start time
Sim.StartTime = SimTime;
Sim.nSimSteps = CLC.Con.horizon_prediction;

% Step 1 Generate the rajectory
[tr, ~] = mr_1_tr(t_n, x, cu, rn, CLC.g0);

if CLC.Con.horizon_prediction > CLC.Con.horizon_action
    Con.YawData = [t', tr;
        SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];
else
    Con.YawData = [t', tr];
end

J = cost_function(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
end












function [T_tmp, Con_tmp, CLC_tmp] = reduceT(T, Con, CLC, group)
% Reduces the wind farm to the turbines that are relevant to the given
% optimization
ind         = zeros(1,T.nT);
ind(group)  = 1:sum(group);
ind_OP      = (repmat(1:T.nOP,sum(group),1) + T.StartI(group)' - 1)';

T_tmp.posBase   = T.posBase(group,:);
T_tmp.nT        = sum(group);
T_tmp.nOP       = T.nOP;
T_tmp.posNac    = T.posNac(group,:);
T_tmp.D         = T.D(group);
T_tmp.Names_OP  = T.Names_OP;
T_tmp.Names_WF  = T.Names_WF;
T_tmp.Names_T   = T.Names_T;
T_tmp.StartI    = (0:T_tmp.nT-1)*T_tmp.nOP + 1;

T_tmp.States_OP = T.States_OP(ind_OP(:),:);
T_tmp.States_WF = T.States_WF(ind_OP(:),:);
T_tmp.States_T  = T.States_T(ind_OP(:),:);

T_tmp.dep       = T.dep(group);
for iD = 1:length(T_tmp.dep)
    T_tmp.dep{iD} = ind(T_tmp.dep{iD});
end

% ======= Con_tmp =======
ind_g = 1:T.nT;
ind_g = ind_g(group);

Con_tmp.Yaw         = Con.Yaw;
Con_tmp.YawData     = Con.YawData(:,[1,ind_g+1]);
Con_tmp.yawRangeMin = CLC.Con.yawRangeMin;
Con_tmp.yawRangeMax = CLC.Con.yawRangeMax;
Con_tmp.tanhYaw     = CLC.Con.tanhYaw;
% ======= CLC_tmp =======
CLC_tmp     = CLC;
CLC_tmp.gd  = CLC.gd(group);
CLC_tmp.g0  = CLC.g0(group);
end