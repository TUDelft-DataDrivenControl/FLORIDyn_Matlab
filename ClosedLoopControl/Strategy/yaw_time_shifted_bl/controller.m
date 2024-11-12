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

%% Step 1: Retrieve the optimisation problem
Uinfty = mean(T.States_WF(T.StartI,1)); % Double check
[T_star,combined,~,ph,solve_order] = ...
    calcT(T, Uinfty, CLC.Con.horizon_action, Sim.TimeStep);

CLC.Con.horizon_prediction = ph;
%% Step 2: Solve optimisation problems
x0_fixed = false(size(CLC.x0));

% Constants
% Time line (+1 to account for current time step)
t   = linspace(SimTime, SimTime + (CLC.Con.horizon_action-1) * Sim.TimeStep,...
    CLC.Con.horizon_action);

if debug; figure(123456); clf; db_counter = 1; end

tic
x0 = CLC.x0;
options = optimoptions('fmincon','MaxIterations',100);
f = @(x) cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,CLC.Set.yawRateLimit,SimTime,t);
mycon =@(x) mycon_wrapper(x,T,Con,CLC,Sim.TimeStep);

% orien = wind dir - yaw;
mean_dir = mean(T.States_WF(T.StartI,2));                           % Needs change for 360 -> 0 transition
ub_orien = mean_dir - CLC.Set.yawRangeMin; % CLC.Set.yawRangeMin<=0
lb_orien = mean_dir - CLC.Set.yawRangeMax; % CLC.Set.yawRangeMax>=0

x = fmincon(f,x0,[],[],[],[],...
    ones(length(x0),1)*lb_orien,...
    ones(length(x0),1)*ub_orien,...
            mycon,options);

CLC.x0 = x;
time_elapsed = toc;
disp(['Naive controller solve time: ' num2str(time_elapsed) ' s for ' num2str(length(x)) ' variables'])

%% Plot random solution space
debug = false;
if debug
    e1 = rand(size(x))-.5;
    e1 = e1./sqrt(sum(e1.^2));
    
    E = null(e1');
    ie2 = randi(size(E,2));
    e2 = E(:,ie2)./sqrt(sum(E(:,ie2).^2));
    
    e_steps = linspace(-30,30,31);
    [E1,E2] = meshgrid(e_steps,e_steps);
    J_E = zeros(size(E1));
    for ie = 1:length(E1(:))
        xe = E1(ie)*e1 + E2(ie)*e2 + x;%ones(size(x))*.5;
        
        if sum(abs(mean_dir - xe)>30)>0
            J_E(ie)=0; 
        else
            J_E(ie) = f(xe);
        end
    end
    J_E(J_E==0)=nan;
    figure
    contourf(E1,E2,J_E,31,'EdgeColor','w')
    hold on
    scatter(0,0,20,"white",'filled');
    hold off
    axis equal
    grid on
    xlabel('e_1')
    ylabel('e_2')
    colormap(viridis(31))

    ax = gca;
    ax.YColor = 'w';
    ax.XColor = 'w';
    ax.Color = 'k';
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
%CLC.x0(Con.shift.relevant_j==1) = x;
CLC.x0 = x;
tr = reshape(CLC.x0,CLC.Con.horizon_action,T.nT);

Con.YawData = [t', tr;
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];

J = cost_function(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
end

function [c,ceq] = mycon_wrapper(x,T,Con,CLC,DeltaT)
ceq = 0;

limyaw = CLC.Set.yawRateLimit * DeltaT;
%CLC.x0(Con.shift.relevant_j==1) = x;
CLC.x0 = x;

tr = [(T.States_WF(T.StartI,2) - T.States_T(T.StartI,2))';...
    reshape(CLC.x0,CLC.Con.horizon_action,T.nT)];
tr_diff = reshape(tr(2:end,:)-tr(1:end-1,:),[],1);

c = abs(tr_diff) - limyaw;
end

%% NOTES
% [ ] fixed and not fixed parameters are ignored at the moment