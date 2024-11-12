function Con = controller(T, Wind, Sim, Con, Vis, paramFLORIDyn, paramFLORIS,CLC, SimTime)
%CONTROLLER determines the control inputs for the near future and beyond.
%   It utilizes a given cost-function to evaluate the 
%% Bezier_yaw
% This controller is based on Bezier-Curves, which allow a smooth
% transition between distinct setpoints. These setpoints are subject to the
% optimization.
% Bezier-curves further allow the user to calculate the derivative over
% time, which allows us to enforce a yaw-rate limit (or axial induction)

%% Task 1 Map out the solution space

% Calculate the rate limit
%   The rate limit for the normed space is calculated by multiplying the
%   rad/s with the time window, which leads to rad/[-].
rate_lim = CLC.Set.yawRateLimit * Sim.TimeStep * CLC.Con.horizon_action;
%   In future versions of the code, one could also change the action time
%   for each turbine, allowing for different action horizons.

%% Generate input for turbines
% test coord \in [0,1]
test_coord = [CLC.test;.5,0.5];
yaw_n = BC_3_2dof(rate_lim,test_coord);

% Convert the normalized yaw angles into a trajectory
[t_o,yaw] = BC_3_2dof_trajectory(...
    ones(T.nT,1)*Sim.TimeStep * CLC.Con.horizon_action, ...
    Sim.TimeStep, yaw_n, T.States_WF(T.StartI,2)); % TODO Yaw_0 is WF & yaw -> combine!

% Write the yaw trajectory into the control yaw array
Con.YawData = [SimTime+t_o, yaw];

%% Debugging
%debug_plot_trajectories(t_o,yaw,nT);
end

function debug_plot_trajectories(t_o,yaw,nT)
figure
for iT= 1:nT
    plot(t_o,yaw(:,iT))
    if iT==1;hold on;end
end
hold off
grid on
end