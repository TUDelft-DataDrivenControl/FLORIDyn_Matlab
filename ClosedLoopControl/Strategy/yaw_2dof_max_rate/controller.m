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
%% Max-rate
% This controller steers the turbines with the maximum yaw speed, or does
% not change their angle.

%% Check if controller should be active
if SimTime<CLC.Time.StartTime
    return 
end

%% Predict Wind conditions
Wind_Con = Wind;


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
options = optimoptions('fmincon','MaxIterations',10,'UseParallel',true);%,'Algorithm','sqp');


f = @(x) cost_function_wrapper(x,T,Wind_Con,Sim,Con,Vis,paramFLORIDyn,...
    paramFLORIS,CLC,rate_lim,SimTime,t_n,t,cu);
tic
x = fmincon(f,x0,[],[],[],[],zeros(CLC.Con.DOFpT*T.nT,1),...
    ones(CLC.Con.DOFpT*T.nT,1),[],options);
disp(toc)

%x = x';
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