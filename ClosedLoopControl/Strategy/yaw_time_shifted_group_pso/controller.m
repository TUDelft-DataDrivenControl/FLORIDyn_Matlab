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
debug = false;


%% Calculate the rate limit
rate_lim = CLC.Set.yawRateLimit * Sim.TimeStep;


%% Step 1: Retrieve the optimisation problem
Uinfty = mean(T.States_WF(T.StartI,1));
[T_star,...
    combined,...
    ~,...
    ph,...
    solve_order,...
    t2o,...
    t2g,...
    opt_ph] = ...
        calcT(T, Uinfty, CLC.Con.horizon_action, Sim.TimeStep);

t2o_time = zeros(size(t2o(:,1)));
t2g_time = zeros(length(solve_order),1);
t2o_time_i = 1;

CLC.Con.horizon_prediction = ph;
%% Step 2: Solve optimisation problems
x0_fixed = false(size(CLC.x0));

% Constants
% Time line (+1 to account for current time step)
t   = linspace( ...
    SimTime, ...
    SimTime + (CLC.Con.horizon_action-1) * Sim.TimeStep,...
    CLC.Con.horizon_action);

turbines = 0:T.nT-1;
% Vector that links each entry of x to the turbine it belongs to
x2t = repmat(turbines,CLC.Con.horizon_action,1); x2t = x2t(:);

% if debug; figure(123456); clf; db_counter = 1; end

for iG = 1:length(solve_order) % <- can be parallelised
    % Each independent optimisation problem 
    for iO = 1:length(solve_order{iG})
        tic;

        subset_x    = combined(solve_order{iG}(iO),:)>0;
        subset_t    = t2g(iG,:)>0;

        % Get subset of x that is optimized in this simulation
        x_opt = CLC.x0(subset_x);
        
        % Get subset of the turbines that are simulated in this simulation
        t_opt = turbines(subset_t);
        
        % Get the prediction horizon of the optimization
        ph_opt = opt_ph(solve_order{iG}(iO));
        
        % Get the vector the reduced simulation output needs to be multiplied
        % by:
        T_star_opt = reshape(T_star(solve_order{iG}(iO),:),ph,[]);
        %   Dim 1: Shorten to prediction horizon
        %   Dim 2: Select only relevant turbines
        T_star_opt = T_star_opt(1:ph_opt, t2g(iG,:)>0);
        T_star_opt = T_star_opt(:)';
        
        % Reduce wind farm to relevant turbines and data
        [T_tmp, Con_tmp, CLC_tmp] = reduceT(T, Con, CLC, t2g(iG,:)>0);
        Con_tmp.shift.relevant_t  = T_star_opt;


        % ====================== Optimize =============================== %
        disp(['Simulating the turbine(s) ' num2str(t_opt) ', ' ...
            'actuating turbine(s) ' num2str(unique(x2t(combined(solve_order{iG}(iO),:)>0))')...
            ', optimizing ' num2str(length(x_opt)) ...
            ' values across ' num2str(opt_ph(solve_order{iG}(iO))) ' time steps'])

        % Create a cost function
        % Build cost function wrapper that incorporates the reduced x into x
        f = @(x) cost_function_wrapper2(x,T_tmp,Wind,Sim,Con_tmp,Vis, ...
            paramFLORIDyn, paramFLORIS, CLC_tmp, ph_opt, ...
            SimTime, subset_x, subset_t, t);


        options = optimoptions('particleswarm',...
            'SwarmSize',100,...
            'MaxIterations',CLC.max_iterations,...
            'UseParallel',true,...
            'Display','iter');
            

        x = particleswarm(f,length(x_opt),...
            -ones(length(x_opt),1)*rate_lim,...
            ones(length(x_opt),1)*rate_lim,options);


        % ===================== Store result ============================ %
        CLC.x0(subset_x)    = x;
        x0_fixed(subset_x)  = true;
        disp([sum(x0_fixed)*100 ' % of the dof optmized.'])

        t2o_time(t2o_time_i) = toc;
        t2g_time(iG)         = t2g_time(iG) + t2o_time(t2o_time_i);

        disp(['Optimisation problem ' num2str(iG) '-' ...
            num2str(solve_order{iG}(iO)) ' took ' num2str(debug_time)...
            ' s to solve ' num2str(debug_nx) ' variables out of ' ...
            num2str(length(CLC.x0)) '.'])

        t2o_time_i           = t2o_time_i + 1;
    end
end

opt_time = sum(t2o_time);
disp(['The optimization took ' num2str(opt_time) 's total.'])




%% Store optimised trajectory in Con.YawData
% Generate the rajectory & store
tr      = reshape(CLC.x0,CLC.Con.horizon_action,T.nT);
tr(1,:) = tr(1,:) + reshape(CLC.g0,1,[]);
tr      = cumsum(tr);

Con.YawData = [t', tr;
    SimTime + CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];

CLC.g0 = Con.YawData(CLC.Time.nS,2:end)';


%% Store data for postprocessing
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

    % [~,M] = cost_function(T,Wind,Sim_tmp,Con,Vis,paramFLORIDyn,paramFLORIS);
    % writetable(M,...
    %     [folder_name filesep pad(num2str(SimTime),7,"left",'0') '_M.csv'])
end


end





function J = cost_function_wrapper2(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC, ph_opt, SimTime, subset_x, subset_t, t)

% Step 1 Set correct start time & duration
Sim.StartTime = SimTime;
Sim.nSimSteps = ph_opt;


% Step 2 Generate the trajectory
%   Store yaw changes in complete x0
CLC.x0(subset_x) = x;
%   create matrix of yaw changes
tr = reshape(CLC.x0,CLC.Con.horizon_action,[]);
%   reduce to turbines that are present in this simulation
tr = tr(:, subset_t);
%   add current orientation as offset (.g0 was reduced with the farm)
tr(1,:) = tr(1,:) + reshape(CLC.g0,1,[]);
tr      = cumsum(tr);
%   store trajectories
Con.YawData = [t', tr; ...
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];



% Step 3 Call the cost function
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

