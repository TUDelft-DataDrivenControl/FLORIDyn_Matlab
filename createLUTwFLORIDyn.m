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

%% Create a Steady state LUT with FLORIDyn


function lut = createLUTwFLORIDyn(Dirs, pathToSimulation)

% Load data from the simulation
% Reset the Matlab Path and load essential paths & the simulation path
addPaths;

% Get the settings for the wind field, visualization, controller and Sim.
[Wind, Vis, Sim, Con] = setup();

%Sim.EndTime = Overwrite_EndTime; % <- Reclict from NAWEA 2023 runs

% Add according functions to the search path
addFLORISPaths;
addFLORIDynPaths;

% Load linked data
turbProp        = turbineArrayProperties();
paramFLORIS     = parameterFLORIS();
paramFLORIDyn   = parameterFLORIDyn();

% Preprocess loaded data
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);
clear turbProp

lut = zeros(length(Dirs),T.nT);

for iDir = 1:length(Dirs)
    lb = Dirs(iDir) - 30;
    ub = Dirs(iDir) + 30;
    tic
    f =@(x) cost_function_wrapper(x, Dirs(iDir), ...
        T, Wind, Sim, Con, paramFLORIDyn, paramFLORIS, Vis);

    % Optimize
    options = optimoptions('particleswarm',...
        'SwarmSize',100,...
        'MaxIterations',40,...
        'UseParallel',true,...
        'Display','iter');

    x = particleswarm(f,T.nT,...
        ones(T.nT,1).*lb,...
        ones(T.nT,1).*ub,options);
    
    disp(['Total opt. time: ' num2str(toc) ' s.'])
    disp(round(x,1))
    lut(iDir,:) = round(x,1);
end

end


function J = cost_function_wrapper(yaw, dir, T, Wind, Sim, Con, ...
    paramFLORIDyn, paramFLORIS, Vis)

Wind.Dir    = [0, dir; 1000000, dir];
Con.YawData = [0, yaw(:)'; 1000000, yaw(:)'];

% ====== Init simulation
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

[T,M,~,~] = FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

J = - sum(M.("Power generated [MW]")(end-T.nT+1:end));
end

