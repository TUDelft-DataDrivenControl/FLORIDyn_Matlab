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

function [Con, CLC] = controller(T, Wind, Sim, Con, Vis, paramFLORIDyn, paramFLORIS,CLC, SimTime)
%CONTROLLER determines the control inputs for the near future and beyond.
%   It utilizes a given cost-function to evaluate the 
%% Shifted controller
% Controller aims for synchronization betweencontrol action and effect at
% the downstream turbine.

%% Check if controller should be active
if SimTime<CLC.Time.StartTime
    return 
end
debug = false;

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
t   = linspace(SimTime, SimTime + (CLC.Con.horizon_action-1) * Sim.TimeStep,...
    CLC.Con.horizon_action);

if debug; figure(123456); clf; db_counter = 1; end

for iO = 1:length(solve_order) % <- can be parallelised
    % Each independent optimisation problem 
    for iiO = 1:length(solve_order{iO})
        tic;

        % Reduce wind farm to relevant turbines and data
        [T_tmp, Con_tmp, CLC_tmp] = reduceT(T, Con, CLC, t2g(iO,:)>0);
        
        tubine_ids = 1:T.nT;
        T_tmp.tubine_ids = tubine_ids(t2g(iO,:)>0);

        % Derive reduced vector of T_star that extracts the required power
        % generated enties of all power generated during the ph_opt steps
        ph_opt      = opt_ph(solve_order{iO}(iiO));
        T_star_opt  = reshape(T_star(solve_order{iO}(iiO),:),ph,[]);
        T_star_opt  = T_star_opt(1:ph_opt, t2g(iO,:)>0);
        Con.shift.relevant_t     = T_star_opt(:)';
        Con_tmp.shift.relevant_t = Con.shift.relevant_t;

        % Find which elements of x are set in this optimization
        Con.shift.relevant_j     = combined(solve_order{iO}(iiO),:);
        Con_tmp.shift.relevant_j = reshape(Con.shift.relevant_j,[],T.nT);
        Con_tmp.shift.relevant_j = Con_tmp.shift.relevant_j(:, T_tmp.tubine_ids);
        Con_tmp.shift.relevant_j = Con_tmp.shift.relevant_j(:);

        Con_tmp.shift.relevant_j_ori = Con.shift.relevant_j;

        % Number of time steps the simulation needs to run
        Sim.nSimSteps = ph_opt;

        % Sim.nSimSteps = ...
        %     find(any(reshape(Con.shift.relevant_t,ph,[]),2),1,'last'); % <- this is how long many steps the simulation needs to run find last 1 per turbine

        %% Set up opt problem
        x0 = CLC.x0(Con.shift.relevant_j==1);
        if debug; debug_nx = length(x0); end
        
        %% Cost function
        % ================== Optimise the trajectories ========================== %
        options = optimoptions('fmincon','Display','iter',...
            'MaxIterations',CLC.max_iterations,...
            'UseParallel',true);
        
        % Creating the cost function wrapper
        f = @(x) cost_function_wrapper(x,T_tmp,Wind,Sim,Con_tmp,Vis,paramFLORIDyn,...
            paramFLORIS,CLC_tmp,CLC.Set.yawRateLimit,SimTime,t);

        % Creating the yawing rate limit wrapper
        mycon =@(x) mycon_wrapper(x,T_tmp,Con_tmp,CLC_tmp,Sim.TimeStep);

        % orien = wind dir - yaw;
        mean_dir = mean(T_tmp.States_WF(T_tmp.StartI,2));                           % Needs change for 360 -> 0 transition 
        ub_orien = mean_dir - CLC.Con.yawRangeMin; % CLC.Set.yawRangeMin<=0
        lb_orien = mean_dir - CLC.Con.yawRangeMax; % CLC.Set.yawRangeMax>=0

        % Run the optimization
        x = fmincon(f,x0,[],[],[],[],...
            ones(sum(Con.shift.relevant_j),1)*lb_orien,...
            ones(sum(Con.shift.relevant_j),1)*ub_orien,...
            mycon,options);
        
        %% Save results
        CLC.x0(Con.shift.relevant_j==1)      = x;
        x0_fixed(Con.shift.relevant_j==1)    = true;
        
        t2o_time(t2o_time_i) = toc;
        t2o_time_i           = t2o_time_i + 1;
        t2g_time(iO)         = t2g_time(iO) + t2o_time(t2o_time_i);

        if debug
            subplot(2,1,1)
            hold on
            imagesc(ones(size(CLC.x0))*db_counter, ...
                1:length(CLC.x0),CLC.x0)
            hold off
            title('Optimised values')
            subplot(2,1,2)
            hold on
            imagesc(ones(size(CLC.x0))*db_counter, ...
                1:length(CLC.x0),x0_fixed)
            hold off
            title('Fixed parameters')

            %% Debug output
            debug_time = toc;
        end
        disp(['Optimisation problem ' num2str(iO) '-' ...
            num2str(solve_order{iO}(iiO)) ' took ' num2str(debug_time)...
            ' s to solve ' num2str(debug_nx) ' variables out of ' ...
            num2str(length(CLC.x0)) '.'])
    end
end

opt_time = sum(t2o_time);
disp(['The optimization took ' num2str(opt_time) 's total.'])

%% Store optimised trajectory in Con.YawData
% Generate the rajectory & store
tr = reshape(CLC.x0,CLC.Con.horizon_action,T.nT);
Con.YawData = [t', tr;
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];

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

%% Shift x0 for each turbine
CLC.x0 = circshift(reshape(CLC.x0,CLC.Con.horizon_action, T.nT), ...
    - CLC.Time.nS);
CLC.x0(end-CLC.Time.nS+1:end,:) = ...
    repmat(CLC.x0(end-CLC.Time.nS,:),CLC.Time.nS,1);
CLC.x0 = CLC.x0(:);
end

function J = cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rn,SimTime,t)

% Correct start time
Sim.StartTime = SimTime;
%Sim.nSimSteps = CLC.Con.horizon_prediction; % <= Already set in outer loop

% Step 1 Generate the trajectory
%[tr, ~] = mr_2_tr(t_n, x, cu, rn, CLC.g0);
CLC.x0(Con.shift.relevant_j==1) = x;
tr = reshape(CLC.x0,CLC.Con.horizon_action,[]);

Con.YawData = [t', tr(:, T.tubine_ids);
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end, T.tubine_ids)];

J = cost_function(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
end

function [c,ceq] = mycon_wrapper(x,T,Con,CLC,DeltaT)
% Ensure that the yawing rate is below the limit
ceq = 0;

limyaw = CLC.Set.yawRateLimit * DeltaT;
CLC.x0(Con.shift.relevant_j_ori==1) = x;

tr_complete = reshape(CLC.x0,CLC.Con.horizon_action,[]);

tr = [(T.States_WF(T.StartI,2) - T.States_T(T.StartI,2))';...
    tr_complete(:,T.tubine_ids)];
tr_diff = reshape(tr(2:end,:)-tr(1:end-1,:),[],1);

c = abs(tr_diff(Con.shift.relevant_j==1)) - limyaw;
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


%% NOTES
% [ ] fixed and not fixed parameters are ignored at the moment