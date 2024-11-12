% ====== Load data from the simulation ======
% Reset the Matlab Path and load essential paths & the simulation path
addPaths;

% Get the settings for the wind field, visualization, controller and Sim.
[Wind, Vis, Sim, Con] = setup();

% ====== Add according functions to the search path
addFLORISPaths;
addFLORIDynPaths;
addEnKFPaths;

% ====== Load linked data
turbProp        = turbineArrayProperties();
paramFLORIS     = parameterFLORIS();
paramFLORIDyn   = parameterFLORIDyn();
CLC             = clc_settings(Sim);
EnKF            = EnKF_settings();
addCLCPaths;

% ====== Preprocess loaded data
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);

% ====== Init simulation
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

% ====== Process and init EnKF Data
[EnKF, Wind] = EnKF_PreProcessing(EnKF, T, Wind, Sim);

% ====== Visualization Settings
Vis_Settings.FlowFieldError.Online  = Vis.FlowField.Error.Online;
Vis.FlowField.Error.Online          = 0;
EnKF_Vis = EnKF_Visualizer(Vis_Settings);

% ====== Init Control parameters ======
CLC.gd = zeros(T.nT,1);
CLC.g0 = T.States_WF(T.StartI,2); % TODO Yaw_0 is WF & yaw

switch CLC.Con.InitMethod
    case "wind dir +- rand"
        CLC.x0 = mean(T.States_WF(T.StartI,2)) + ...
            (2*rand(CLC.Con.DOFpT*T.nT,1)-1);
    case "0 + rand"
        CLC.x0 = rand(CLC.Con.DOFpT*T.nT,1);
    otherwise
        CLC.x0 = rand(CLC.Con.DOFpT*T.nT,1);
end


% ====== Init Output =====
if or(CLC.Con.StoreData, EnKF.StoreData)
    % Case name is date and time
    CLC.CaseName = ['Case' convertStringsToChars(string(datetime(...
        "now",'Format','yyyy_MM_dd_HH_mm_ss')))];
    mkdir([Sim.PathToSim 'Results' filesep CLC.CaseName])
    
    % Automatically move settings to case folder to reproduce
    mkdir([Sim.PathToSim 'Results' filesep CLC.CaseName filesep...
        'Settings'])
    files2copy = {'setup.mlx', ...
        'EnKF_settings.m', ...
        'clc_settings.mlx', ...
        'parameterFLORIDyn.m', 'parameterFLORIS.m'};
    
    for iF = 1:length(files2copy)
        copyfile(...
            [Sim.PathToSim files2copy{iF}],...
            [Sim.PathToSim 'Results' filesep CLC.CaseName ...
            filesep 'Settings' filesep files2copy{iF}])
    end
    clear files2copy iF
end

if isempty(gcp('nocreate'))
    parpool([8,80])
end
clear Vis_Settings turbProp 