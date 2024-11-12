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

% ======================================================================= %
% ================== Optimise the trajectories ========================== %
%options = optimoptions('fmincon','Display','iter','MaxIterations',10);%,'Algorithm','sqp');
options = optimoptions('fmincon','MaxIterations',10);%,'Algorithm','sqp');


f = @(x) cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rate_lim,SimTime,t_n,t,cu);
%% Particle swarm optimisation
options = optimoptions('particleswarm',...
    'SwarmSize',100,...
    'MaxIterations',10,...
    'UseParallel',true,...
    'Display','iter');
x = particleswarm(f,length(x0),...
    zeros(CLC.Con.DOFpT*T.nT,1),...
    ones(CLC.Con.DOFpT*T.nT,1),options);
x = x';
% ======================================================================= %
% MISSING: yaw rate enforcement
% ======================================================================= %

%% Plot random solution space
debug = false;
if debug
    e1 = rand(length(x),1)-0.5;
    e1 = e1./sqrt(sum(e1.^2));
    
    E = null(e1');
    e2 = E(:,randi(size(E,2)));
    e2 = e2./sqrt(sum(e2.^2));
    
    e_steps = linspace(-.5,.5,11);
    [E1,E2] = meshgrid(e_steps,e_steps);
    J_E = zeros(size(E1));
    for ie = 1:length(E1(:))
        xe = E1(ie)*e1 + E2(ie)*e2 + x;%ones(size(x))*.5;
        
        if sum(or(xe>1,xe<0))>0
            J_E(ie)=0; 
        else
            J_E(ie) = f(xe);
        end
    end
    J_E(J_E==0)=nan;
    figure
    contourf(E1,E2,J_E,100,'EdgeColor','white')
    hold on
    scatter(0,0,20,"white",'filled');
    hold off
    axis equal
    grid on
    xlabel('e_1')
    ylabel('e_2')
    colormap(viridis(100))
end

%% Store optimised trajectory in Con.YawData
% Step 1 Generate the rajectory
[tr, CLC.g0] = mr_2_tr(t_n, x, cu, rate_lim, CLC.g0);
CLC.x0 = x;
if CLC.Con.horizon_prediction > CLC.Con.horizon_action
    Con.YawData = [t', tr;
        SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];
else
    Con.YawData = [t', tr];
end
end

function J = cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rn,SimTime,t_n,t,cu)

% Correct start time
Sim.StartTime = SimTime;
Sim.nSimSteps = CLC.Con.horizon_prediction;

% Step 1 Generate the rajectory
[tr, ~] = mr_2_tr(t_n, x, cu, rn, CLC.g0);

if CLC.Con.horizon_prediction > CLC.Con.horizon_action
    Con.YawData = [t', tr;
        SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];
else
    Con.YawData = [t', tr];
end

J = cost_function(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
end