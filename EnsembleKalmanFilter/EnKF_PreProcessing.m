function [EnKF, Wind] = EnKF_PreProcessing(EnKF, T, Wind, Sim)
%ENKF_PREPROCESSING calculates variables based on the EnKF settings and
%preallocates the matrices for the states.

%% Simulation Sections
EnKF.Sim.Sections   = max(floor(Sim.nSimSteps / EnKF.nS),1);
EnKF.Sim.SecDur     = Sim.TimeStep * EnKF.nS;
EnKF.Sim.StartTime  = Sim.StartTime;

%% Limits
if ~isfield(EnKF.Vel,'Limits')
    EnKF.Vel.Limits = false;
end


%% Ensemble States
if EnKF.Vel.Correct
    % Velocity states
    EnKF.States.Vel   = repmat(T.States_WF(:,1),1,EnKF.nE);
    
    % Velocity measurement noise (NOT process noise!)
    if isfield(EnKF.Vel,'MeasurementNoise')
        [EnKF.Vel.C_ee_Vel, EnKF.Vel.C_ee_Vel_Chol] = readCovMatrix(...
            EnKF.Vel.MeasurementNoise,T.nT,'WindVel');
    else
        %EnKF.Vel.C_ee_Vel = readmatrix('WindVelCovariance.csv');
        [EnKF.Vel.C_ee_Vel, EnKF.Vel.C_ee_Vel_Chol] = readCovMatrix(...
            readmatrix('WindVelCovariance.csv'),T.nT,'WindVel');
    end
    for iE = 1:EnKF.nE
        randomError = (randn(T.nOP,T.nT)*EnKF.Vel.C_ee_Vel_Chol);
        EnKF.States.Vel(:,iE) = EnKF.States.Vel(:,iE) + randomError (:);
    end
    
    % Power generated noise
    [EnKF.Output.C_ee_Pow, EnKF.Output.C_ee_Pow_Chol] = readCovMatrix(...
        EnKF.Output.PowNoiseVar,T.nT,'WindVel');

    [EnKF.Vel.C_ee_PowR, EnKF.Vel.C_ee_PowR_Chol] = readCovMatrix(...
        readmatrix('WindPowRootCovariance.csv'),T.nT,'WindDir');
    
end
if EnKF.Dir.Correct
    EnKF.States.Dir = repmat(T.States_WF(:,2),1,EnKF.nE);

    if isfield(EnKF.Dir,'MeasurementNoise')
        [EnKF.Dir.C_ee_Dir, EnKF.Dir.C_ee_Dir_Chol] = readCovMatrix(...
            EnKF.Dir.MeasurementNoise,T.nT,'WindDir');
    else
        [EnKF.Dir.C_ee_Dir, EnKF.Dir.C_ee_Dir_Chol] = readCovMatrix(...
            readmatrix('WindDirCovariance.csv'),T.nT,'WindDir');    
    end
    
    % This line overwrote the Vel.C_ee_Pow settings ... 
    % [EnKF.Vel.C_ee_Pow, EnKF.Vel.C_ee_Pow_Chol] = readCovMatrix(...
    %     readmatrix('WindDirCovariance.csv'),T.nT,'WindDir');
    
    for iE = 1:EnKF.nE
        randomError = (randn(T.nOP,T.nT)*EnKF.Dir.C_ee_Dir_Chol);
        EnKF.States.Dir(:,iE) = EnKF.States.Dir(:,iE) + randomError (:);
    end
end
if EnKF.TI.Correct
    EnKF.States.TI  = repmat(T.States_WF(:,3),1,EnKF.nE);
    [EnKF.TI.C_ee_TI, EnKF.TI.C_ee_TI_Chol] = readCovMatrix(...
        readmatrix('WindTICovariance.csv'),T.nT,'WindTI');
    for iE = 1:EnKF.nE
        randomError = (randn(T.nOP,T.nT)*EnKF.TI.C_ee_TI_Chol);
        EnKF.States.TI(:,iE) = EnKF.States.TI(:,iE) + randomError (:);
    end
end

EnKF.States_OP = repmat(T.States_OP,1,EnKF.nE);
EnKF.nStatesOP = length(T.Names_OP);
EnKF.States_T  = repmat(T.States_T,1,EnKF.nE);
EnKF.nStatesT  = length(T.Names_T);
%% Measurements
EnKF.M = cell(EnKF.nE,1);
EnKF.Interaction = cell(EnKF.nE,T.nT);
EnKF.InteractionNames = 'T_ID, iOP1, wOP1, iOP2, wOP2, w';
EnKF.Output.Pow = zeros(T.nT,EnKF.nE);
%% Flow field inputs
% Allocate fields
EnKF.Wind = Wind;

if EnKF.Vel.Correct
    switch lower(EnKF.Vel.PredicionModel)
        case 'rw'
            Wind = rmfield(Wind,'Vel');
            Wind.Input_Vel = 'EnKF_RW';
            Wind.Vel.CholSig = EnKF.Vel.C_ee_Vel_Chol;
        case 'zoh'
            Wind.Input_Vel = 'EnKF_ZOH';
            Wind.Vel = -1;
        case 'interpolate'
            Wind.Input_Vel = 'EnKF_InterpTurbine';
            Wind.Vel = [[...
                EnKF.Sim.StartTime:Sim.TimeStep:EnKF.Sim.StartTime + ...
                EnKF.Sim.SecDur]',zeros(EnKF.nS+1,T.nT)];
        otherwise
            Wind = rmfield(Wind,'Vel');
            Wind.Input_Vel = 'EnKF_RW';
            Wind.Vel.CholSig = EnKF.Vel.C_ee_Vel_Chol;
    end
    % Enable process noise
    %   Necessary to model the increasing uncertainty in areas where the
    %   values can not be corrected and to encourange correction in areas
    %   where values can be corrected.
    Wind.Pertubation.Vel = true;
    Wind.Pertubation.VelSigma = EnKF.Vel.ProcessNoise;
else
    % Disable process noise
    %   If process noise is on for a state which is not corrected the
    %   values can develop into unreasonable magitudes (negative wind
    %   speeds etc.)
    Wind.Pertubation.Vel = false;
end


if EnKF.Dir.Correct
    switch lower(EnKF.Dir.PredicionModel)
        case 'rw'
            Wind = rmfield(Wind,'Dir');
            Wind.Input_Dir = 'EnKF_RW';
            Wind.Dir.CholSig = EnKF.Dir.C_ee_Dir_Chol;
        case 'zoh'
            Wind = rmfield(Wind,'Dir');
            Wind.Input_Dir = 'EnKF_ZOH';
            Wind.Dir = -1;
        case 'interpolate'
            Wind.Input_Dir = 'EnKF_InterpTurbine';
            Wind.Dir = [[...
                EnKF.Sim.StartTime:Sim.TimeStep:EnKF.Sim.StartTime + ...
                EnKF.Sim.SecDur]',zeros(EnKF.nS+1,T.nT)];
        otherwise
            Wind = rmfield(Wind,'Dir');
            Wind.Input_Dir = 'EnKF_RW';
            Wind.Dir.CholSig = EnKF.Dir.C_ee_Dir_Chol;
    end
    % Enable process noise
    %   Necessary to model the increasing uncertainty in areas where the
    %   values can not be corrected and to encourange correction in areas
    %   where values can be corrected.
    Wind.Pertubation.Dir = true;
    Wind.Pertubation.DirSigma = EnKF.Dir.ProcessNoise;
else
    % Disable process noise
    %   If process noise is on for a state which is not corrected the
    %   values can develop into unreasonable magitudes (negative wind
    %   speeds etc.)
    Wind.Pertubation.Dir = false;
end


if EnKF.TI.Correct
    switch lower(EnKF.TI.PredicionModel)
        case 'rw'
            Wind = rmfield(Wind,'TI');
            Wind.Input_TI = 'EnKF_RW';
            Wind.Dir.CholSig = EnKF.TI.C_ee_TI_Chol;
        case 'zoh'
            Wind.Input_TI = 'EnKF_ZOH';
            Wind.Vel = -1;
        case 'interpolate'
            Wind.Input_TI = 'EnKF_InterpTurbine';
            Wind.TI = [[...
                EnKF.Sim.StartTime:Sim.TimeStep:EnKF.Sim.StartTime + ...
                EnKF.Sim.SecDur]',zeros(EnKF.nS+1,T.nT)];
        otherwise
            Wind = rmfield(Wind,'TI');
            Wind.Input_TI = 'EnKF_RW';
            Wind.Dir.CholSig = EnKF.TI.C_ee_TI_Chol;
    end
    % Enable process noise
    %   Necessary to model the increasing uncertainty in areas where the
    %   values can not be corrected and to encourange correction in areas
    %   where values can be corrected.
    Wind.Pertubation.TI = true;
    Wind.Pertubation.TISigma = EnKF.TI.ProcessNoise;
else
    % Disable process noise
    %   If process noise is on for a state which is not corrected the
    %   values can develop into unreasonable magitudes (negative wind
    %   speeds etc.)
    Wind.Pertubation.TI = false;
end
end