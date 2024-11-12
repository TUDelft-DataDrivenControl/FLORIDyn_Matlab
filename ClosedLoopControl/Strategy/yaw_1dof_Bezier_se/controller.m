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

% Constants
% Time line (+1 to account for current time step)
t_n = linspace(0,1,CLC.Con.horizon_action+1);
t   = linspace(SimTime, SimTime + CLC.Con.horizon_action * Sim.TimeStep,...
    CLC.Con.horizon_action+1);

% Control Update ratio
cu = CLC.Time.nS/CLC.Con.horizon_action;

% ======================================================================= %
% ================== Optimise the trajectories ========================== %
%options = optimoptions('fmincon','Display','iter','MaxIterations',10);%,'Algorithm','sqp');
options = optimoptions('fmincon','MaxIterations',10);%,'Algorithm','sqp');

f = @(x) cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rate_lim,SimTime,t_n,t,cu);
x = fmincon(f,x0,[],[],[],[],zeros(T.nT,1),ones(T.nT,1),[],options);
% ======================================================================= %
% MISSING: yaw rate enforcement
% ======================================================================= %

%% Store optimised trajectory in Con.YawData
% Step 1 Generate the rajectory
[tr, CLC.g0, CLC.gd] = BC_3_1dof_trajectory(t_n,x,CLC.gd,rate_lim,cu,CLC.g0);
CLC.x0 = x;
Con.YawData = [t', tr;
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];
end

function J = cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rn,SimTime,t_n,t,cu)

% Correct start time
Sim.StartTime = SimTime;
Sim.nSimSteps = CLC.Con.horizon_prediction;

% Step 1 Generate the rajectory
[tr, ~, ~] = BC_3_1dof_trajectory(t_n,x,CLC.gd,rn,cu,CLC.g0);

Con.YawData = [t', tr;
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:)];

J = cost_function(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
end