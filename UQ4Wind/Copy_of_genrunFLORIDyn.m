% This file is used to run and initialize all settings in FLORIDyn in order 
% to perform both sensitivity analysis and Bayesian calibration

% Output Y of FLORIDyn for input X
% X comprises of parameters defined in uq_sensitivity.m and uq_calibration
function Y = genrunFLORIDyn(X)
cd ..
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
paramFLORIS.d       = X(9);
paramFLORIS.eta     = X(10);
paramFLORIS.p_p     = X(11);

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
[T,M,Vis] = FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
disp('done')

%%



main_addPaths;

%% Set controller type
% Setting for the contoller
% Control.Type:
%   'SOWFA_greedy_yaw'  -> Uses SOWFA yaw angles and a greedy controller
%                           for C_T and C_P based on lookup tables and the
%                           wind speed (needs additional files)
%   'SOWFA_bpa_tsr_yaw' -> Uses SOWFA yaw angles, blade-pitch-angles and
%                           tip-speed-ratio (needs additional files)
%   'FLORIDyn_greedy'   -> A greedy controller based on lookup tables and
%                           the wind speed (no additional files)
%
% Control.init:
%   Set to true if you are starting a new simulation, if you are copying
%   the states from a previous simulation, set to false.

Control.Type = 'SOWFA_greedy_yaw';
Control.init = true;

%% User Choice

%global n;
n = 3;                                        % 3T case
filename = 'threeDTU10MW_';

scenario = 'const';                           %   Currently implemented scenarios:
                                              %       'const'                     -> Constant wind speed, direction and
                                              %                                       amb. turbulence
                                              %       '+60DegChange'              -> 60 degree wind angle change after
                                              %                                       300s (all places at the same time)


                                              %       'Propagating40DegChange'    -> Propagating 40 degree wind angle
                                              %                                       change starting after 300s

wind_angle = 0;
simulation_duration = 1000; 
freespeed = true;
interaction = true;
posterior_measure_factor = 2000;
alpha_z = 0.08;
wind_speed = 8;
ambient_turbulence = 0.062;

%% Set path to SOWFA files
% To run this, modify the SOWFA output files to have the ending .csv they
% are expected to be avaiable under i.e. [file2val 'generatorPower.csv']
%
% Two control options are implemented:
%   1) Run greedy control
%       Needed files:
%       'nacelleYaw.csv'
%   2) Calculate Ct and Cp based on blade pitch and tip speed ratio
%       Needed files:
%       'nacelleYaw.csv','generatorPower.csv',
%       'bladePitch.csv','rotorSpeedFiltered.csv'
%       ATTENTION!
%       The SOWFA file 'bladePitch.csv' has to be modified to say
%           0     instead of     3{0}
%       Search & delete all "3{" and "}"
%
% Needed for plotting:
% 'generatorPower.csv'

file2val = filename; % 3T case
LoadSOWFAData;

Sim.X = X;

%% Load Layout
%   Load the turbine configuration (position, diameter, hub height,...) the
%   power constants (Efficiency, p_p), data to connect wind speed and
%   power / thrust coefficient and the configuration of the OP-chains:
%   relative position, weights, lengths etc.
%
%   Currently implemented Layouts
%       'oneDTU10MW'    -> one turbine
%       'twoDTU10MW'    -> two turbines at 900m distance
%       'nineDTU10MW'   -> nine turbines in a 3x3 grid, 900m dist.
%       'threeDTU10MW'  -> three turbines in 1x3 grid, 5D distance
%       'fourDTU10MW'   -> 2x2 grid
%
%   Chain length & the number of chains can be set as extra vars, see
%   comments in the function for additional info.

[T,fieldLims,Pow,VCpCt,chain] = loadLayout(filename,Sim);

%% Load the environment
%   U provides info about the wind: Speed(s), direction(s), changes.
%   I does the same, but for the ambient turbulence, UF hosts constant
%   used for the wind field interpolation, the air density, atmospheric
%   stability etc. The Sim struct holds info about the simulation: Duration
%   time step, various settings. See comments in the function for
%   additional info.
%
%   Currently implemented scenarios:
%       'const'                     -> Constant wind speed, direction and
%                                       amb. turbulence
%       '+60DegChange'              -> 60 degree wind angle change after
%                                       300s (all places at the same time)


%       'Propagating40DegChange'    -> Propagating 40 degree wind angle
%                                       change starting after 300s
%
%   Numerous settings can be set via additional arguments, see the comments
%   for more info.

[U, I, UF, Sim] = loadWindField(scenario,... 
    'windAngle',wind_angle,...
    'SimDuration',simulation_duration,...
    'FreeSpeed',freespeed,...
    'Interaction',interaction,...
    'posMeasFactor',posterior_measure_factor,...
    'alpha_z',alpha_z,...
    'windSpeed',wind_speed,...
    'ambTurbulence',ambient_turbulence);

Sim.X = X;

%% Visulization
% Set to true or false
%   .online:      Scattered OPs in the wake with quiver wind field plot
%   .Snapshots:   Saves the Scattered OP plots, requires online to be true
%   .FlowField:   Plots the flow field at the end of the simulation
%   .PowerOutput: Plots the generated power at the end of the simulation
%   .Console:     Online simulation progress with duration estimation
%                 (very lightweight, does not require online to be true)

Vis.online      = false;
Vis.Snapshots   = false;
Vis.FlowField   = false;
Vis.PowerOutput = false;
Vis.Console     = true;

%% Create starting OPs and build opList
%   Creates the observation point struct (OP) and extends the chain struct.
%   Here, the distribution of the OPs in the wake is set. Currently the
%   avaiable distributions are:
%   'sunflower'         : Recommended distibution with equal spread of the
%                           OPs across the rotor plane.
%   '2D_horizontal'     : OPs in two horizontal planes, silightly above and
%                           below hub height
%   '2D_vertical'       : OPs in two vertical planes, right and left of the
%                           narcelle.

[OP, chain] = assembleOPList(chain,T,'sunflower');

%% Running FLORIDyn

[powerHist,OP,T,chain]=...
    FLORIDyn(T,OP,U,I,UF,Sim,fieldLims,Pow,VCpCt,chain,Vis,Control);

%% Load interpolated SOWFA data

powSOWFA_WPS = importGenPowerFile([file2val 'generatorPower.csv']);

for i = 1:n
    p(:,i)=interp1(...
      powSOWFA_WPS(i:n:end,2) - powSOWFA_WPS(i,2),...
      powSOWFA_WPS(i:n:end,3)/UF.airDen,Sim.TimeSteps);
    
    figure(i);
    plot(powerHist(:,1),powerHist(:,i+1));
    hold on
    
    plot(powerHist(:,1),p(:,i));
    hold on
    
    power_vectorized(:,i) = powerHist(:,i+1);
    %Y = reshape(power_vectorized,1,i*753);
    
end

%Y = [powerHist(:,2);powerHist(:,3);powerHist(:,4)]';       %calculation results

[dimension1,dimension2] = size(power_vectorized);
Y = reshape(power_vectorized,1,dimension1*dimension2);

end
