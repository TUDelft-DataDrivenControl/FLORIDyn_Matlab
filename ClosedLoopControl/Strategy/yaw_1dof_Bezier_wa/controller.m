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
%This controller adapts the final yaw rate of change to be equivalent to
%the rate of change of the wind.
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

% Wind angle rate of change at the end of the time horizon
wa = -(getWindDirT(Wind.Dir,1:T.nT,...
    SimTime + CLC.Con.horizon_action * Sim.TimeStep - Sim.TimeStep) - ...
    getWindDirT(Wind.Dir,1:T.nT,...
    SimTime + CLC.Con.horizon_action * Sim.TimeStep + Sim.TimeStep))/...
    (2 * Sim.TimeStep) * Sim.TimeStep * CLC.Con.horizon_action;

wa = max(-rate_lim,min(rate_lim,wa));

% ======================================================================= %
% ================== Optimise the trajectories ========================== %
%options = optimoptions('fmincon','Display','iter','MaxIterations',10);%,'Algorithm','sqp');
%options = optimoptions('fmincon','Display','iter','MaxIterations',10);%,...
    %'PlotFcn',{'optimplotx','optimplotfval','optimplotstepsize'});
options = optimoptions('fmincon','MaxIterations',10);%,'Algorithm','sqp');
f = @(x) cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rate_lim,SimTime,t_n,t,cu,wa);


% ///////////// Brute force the cost //////////////////
% x_t = x0;
% t_var = linspace(0,1,11);
% for ix = 1:11
%     x_t(1) = t_var(ix);
%     for iy = 1:11
%         x_t(2) = t_var(iy);
%         CLC.costFunc(iy,ix,iCLC) = f(x_t);
%     end
% end
% /////////////////////////////////////////////////////

x = fmincon(f,x0,[],[],[],[],zeros(T.nT,1),ones(T.nT,1),[],options);

% % ///// store optimized x
% CLC.final_x(1:2,iCLC) = x;

% === Test plot for cost and the chosen value
% figure(23);
% contourf(CLC.costFunc(:,:,iCLC));
% hold on
% scatter(x(1)*10+1,x(2)*10+1,'r','filled');
% scatter(x0(1)*10+1,x0(2)*10+1,'r');
% hold off
% pause(0.1)
% ======================================================================= %
% MISSING: yaw rate enforcement
% ======================================================================= %

%% Store optimised trajectory in Con.YawData
% Step 1 Generate the rajectory
[tr, CLC.g0, CLC.gd] = BC_3_1dof_trajectory(t_n,x,CLC.gd,rate_lim,...
    cu,wa,CLC.g0);
CLC.x0 = x;
Con.YawData = [t', tr;
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:) + ...
    wa/(Sim.TimeStep * CLC.Con.horizon_action) *...
    (CLC.Con.horizon_prediction-CLC.Con.horizon_action)*Sim.TimeStep];
end

function J = cost_function_wrapper(x,T,Wind,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rn,SimTime,t_n,t,cu,wa)

% Correct start time
Sim.StartTime = SimTime;
Sim.nSimSteps = CLC.Con.horizon_prediction;

% Step 1 Generate the rajectory
[tr, ~, ~] = BC_3_1dof_trajectory(t_n,x,CLC.gd,rn,cu,wa,CLC.g0);

Con.YawData = [t', tr;
    SimTime+CLC.Con.horizon_prediction*Sim.TimeStep, tr(end,:) + ...
    wa/(Sim.TimeStep * CLC.Con.horizon_action) *...
    (CLC.Con.horizon_prediction-CLC.Con.horizon_action)*Sim.TimeStep];

J = cost_function(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
end