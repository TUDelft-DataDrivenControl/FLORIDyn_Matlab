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
debug = true;

%% Step 1: Retrieve the optimisation problem
Uinfty = mean(T.States_WF(T.StartI,1));
[T_star,combined,~,ph,solve_order,~,~] = ...
    calcT(T, Uinfty, CLC.Con.horizon_action, Sim.TimeStep);

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
        if debug; tic; end
        % Dependent opt problems, solved in order
        Con.shift.relevant_t = T_star(solve_order{iO}(iiO),:);
        Con.shift.relevant_j = combined(solve_order{iO}(iiO),:);
        % OR combine 

        Sim.nSimSteps = ...
            find(any(reshape(Con.shift.relevant_t,ph,[]),2),1,'last'); % <- this is how long many steps the simulation needs to run find last 1 per turbine

        %% Set up opt problem
        x0 = CLC.x0(Con.shift.relevant_j==1);
        if debug; debug_nx = length(x0); end
        

        %% Cost function
        % ================== Optimise the trajectories ========================== %
        %options = optimoptions('fmincon','Display','iter','MaxIterations',10);%,'Algorithm','sqp');
        options = optimoptions('fmincon','MaxIterations',100);%,'Algorithm','sqp');
        
        f = @(x) cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
            paramFLORIS,CLC,CLC.Set.yawRateLimit,SimTime,t);
        mycon =@(x) mycon_wrapper(x,T,Con,CLC,Sim.TimeStep);

        % orien = wind dir - yaw;
        mean_dir = mean(T.States_WF(T.StartI,2));                           % Needs change for 360 -> 0 transition 
        ub_orien = mean_dir - CLC.Set.yawRangeMin; % CLC.Set.yawRangeMin<=0
        lb_orien = mean_dir - CLC.Set.yawRangeMax; % CLC.Set.yawRangeMax>=0

        x = fmincon(f,x0,[],[],[],[],...
            ones(sum(Con.shift.relevant_j),1)*lb_orien,...
            ones(sum(Con.shift.relevant_j),1)*ub_orien,...
            mycon,options);
        
        %% Save results
        CLC.x0(Con.shift.relevant_j==1)      = x;
        x0_fixed(Con.shift.relevant_j==1)    = true;

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
            disp(['Optimisation problem ' num2str(iO) '-' ...
                num2str(solve_order{iO}(iiO)) ' took ' num2str(debug_time)...
                ' s to solve ' num2str(debug_nx) 'variables out of ' ...
                num2str(length(CLC.x0)) '.'])
        end
    end
end

%% Store optimised trajectory in Con.YawData
% Generate the rajectory & store
tr = reshape(CLC.x0,CLC.Con.horizon_action,T.nT);
Con.YawData = [t', tr;
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];
end

function J = cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rn,SimTime,t)

% Correct start time
Sim.StartTime = SimTime;
%Sim.nSimSteps = CLC.Con.horizon_prediction; % <= Already set in outer loop

% Step 1 Generate the rajectory
%[tr, ~] = mr_2_tr(t_n, x, cu, rn, CLC.g0);
CLC.x0(Con.shift.relevant_j==1) = x;
tr = reshape(CLC.x0,CLC.Con.horizon_action,T.nT);

Con.YawData = [t', tr;
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];

J = cost_function(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
end

function [c,ceq] = mycon_wrapper(x,T,Con,CLC,DeltaT)
ceq = 0;

limyaw = CLC.Set.yawRateLimit * DeltaT;
CLC.x0(Con.shift.relevant_j==1) = x;

tr = [(T.States_WF(T.StartI,2) - T.States_T(T.StartI,2))';...
    reshape(CLC.x0,CLC.Con.horizon_action,T.nT)];
tr_diff = reshape(tr(2:end,:)-tr(1:end-1,:),[],1);

c = abs(tr_diff(Con.shift.relevant_j==1)) - limyaw;
end

%% NOTES
% [ ] fixed and not fixed parameters are ignored at the moment