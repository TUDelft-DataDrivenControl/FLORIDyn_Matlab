%% Add paths
% Reset paths
restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED

%disp(['Current folder:' pwd])
% Add basic paths
addpath(genpath(['.' filesep 'Visualization']));
addpath(['.' filesep 'FLORIDynCL']);
addpath(['.' filesep 'Data' filesep 'TurbineData']);
addpath(['.' filesep 'Data' filesep 'StoreData']);
addpath(['.' filesep 'Correction' filesep 'Functions'])
addpath(['.' filesep 'EnsembleKalmanFilter'])
addpath(['.' filesep 'ClosedLoopControl'])
addpath(['.' filesep 'Correction' filesep 'GetData'])

% Load simulation folder
addpath(genpath(['.' filesep 'Simulations' filesep pathToSimulation]));
rmpath(genpath(['.' filesep 'Simulations' filesep pathToSimulation filesep 'Results']))