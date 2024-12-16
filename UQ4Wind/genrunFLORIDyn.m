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

% This file is used to run and initialize all settings in FLORIDyn in order 
% to perform both sensitivity analysis and Bayesian calibration

% Output Y of FLORIDyn for input X
% X comprises of parameters defined in uq_sensitivity.m and uq_calibration
function Y = genrunFLORIDyn(X)
%% MainFLORIDyn Center-Line model
% Improved FLORIDyn approach over the gaussian FLORIDyn model

%% Load data from the simulation
% Get the settings for the wind field, visualization, controller and Sim.
[Wind, Vis, Sim, Con] = setup();

% Add according functions to the search path
addFLORISPaths;

% Load linked data
turbProp        = turbineArrayProperties();
paramFLORIDyn   = parameterFLORIDyn();

%% Assign X
%paramFLORIS     = parameterFLORIS();
paramFLORIS.alpha   = X(1);
paramFLORIS.beta    = X(2);
paramFLORIS.k_a     = X(3);
paramFLORIS.k_b     = X(4);
paramFLORIS.k_fa    = X(5);
paramFLORIS.k_fb    = X(6);
paramFLORIS.k_fc    = X(7);
paramFLORIS.k_fd    = X(8);
paramFLORIS.eta     = X(10);
paramFLORIS.p_p     = X(11);
paramFLORIS.airDen  = 1.225;
paramFLORIDyn.advectionFactor = X(9);

%%
global m sim_p
if mod(sim_p,5)==0
    Vis.CL.Store = true;
    Vis.CL.Run = num2str(sim_p);
else
    Vis.CL.Store = false;
end
%% Preprocess loaded data
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);
clear turbProp

%% ====== Init simulation
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

%% ============ RUN SIMULATION ============
% T := Simulation state (OP states, Turbine states, wind field states(OPs))
% M := Measurements from the simulation (T_red, T_addedTI, T_Ueff, T_pow)
tic
[~,M,~] = FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
tmeas = toc;

sim_p = sim_p + 1;
message = [...
    'Progess: ' num2str(round(100*sim_p/m)) ' %%\n'...
    'Remaining time: ' timeEst(tmeas,m-sim_p) ' (hh:mm:ss)'];
dispstat(sprintf(message)) % Can be deleted

%% Plot data
% Change boolean as needed
plotPower = false;
plotUfree = true;
offset = 1000;

timeFDyn = M.("Time [s]") - M.("Time [s]")(1);
if plotPower
    plotQuantityOfInterest(timeFDyn,...
        M.("Power generated [MW]"),T.nT,...
        110,'Power generated [MW]',sim_p/m) 
    Y = reshape(M.("Power generated [MW]")(M.("Time [s]")>=offset),1,[]);
elseif plotUfree
    plotQuantityOfInterest(timeFDyn,...
        M.("Free Wind Speed [ms^-1]"),T.nT,...
        111,'Free Wind Speed [ms^-1]',sim_p/m)
    Y = reshape(M.("Free Wind Speed [ms^-1]"),1,[]);
else
    error('genrunFLORIDyn: No output selected!')
end

end
