function WSE = WSEParameters(nT,pathSimulation,deltaT)
%% Settings for the I & I wind speed estimator
try
    rotorSpeed  = ...
        importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_rotorSpeed.csv']);
catch
    rotorSpeed  = ...
        importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_rotorSpeedFiltered.csv']);
end
bladePitch  = importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_bladePitch.csv']);
genTorque   = importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_generatorTorque.csv']);
nacelleYaw  = importSOWFAFile([pathSimulation 'Data' filesep 'SOWFA_nacelleYaw.csv']);

% Set parameters
WSE.V       = ones(nT,1)*8; % Initial wind speed
WSE.Vinit   = WSE.V;
WSE.gamma   = 5; %20; %5
WSE.beta    = 0; %40; %0
WSE.omega   = rotorSpeed((1:nT),3)* pi/30;
WSE.Ee      = ones(nT,1)*0;

% Time during which the Initial wind speed will be used while the Estimator
% converges
WSE.Offset  = 100; % [s]

% Save turbine and Sim properties
WSE.T_prop  = estimator_dtu10mw();
WSE.dt_SOWF = rotorSpeed(nT+1,2)-rotorSpeed(1,2);
WSE.dt_FDyn = deltaT;
WSE.nT      = nT;

% Save data
WSE.rotorSpeed = rotorSpeed;
WSE.bladePitch = bladePitch;
WSE.genTorque  = genTorque;
WSE.nacelleYaw = nacelleYaw;

end