function [T, Wind, Sim, Con, paramFLORIS, Vis] = prepareSimulation(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim, Vis)
%PREPARESIMULATION Loads the data from the linked folders and generates
%intial conditions

loadDataWarnings = {};
%% Prepare wind
% Load data necessary to run the wind field
% ============= Velocity ============= 
switch Wind.Input_Vel
    case 'I_and_I'
        Wind.Vel.WSE = WSEParameters(...
            size(turbProp.Pos,1),Sim.PathToSim,Sim.TimeStep);
        Wind.Vel.TimePrev  = Sim.StartTime;
        Wind.Vel.StartTime = Sim.StartTime;
    case 'Interpolation'
        try
            Wind.Vel = readmatrix('WindVel.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'WindVel.csv',2,[],...
                [Sim.StartTime,8],[Sim.EndTime,10])
            loadDataWarnings = [loadDataWarnings{:},...
                {['WindVel.csv not found, default created. ' ...
                'Units: [s, ms^-1]']}];
        end
    case 'InterpTurbine'
        try
            Wind.Vel = readmatrix('WindVelTurbine.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'WindVelTurbine.csv',3,size(turbProp.Pos,1),...
                [Sim.StartTime,8],[Sim.EndTime,10])
            loadDataWarnings = [loadDataWarnings{:},...
                {['WindVelTurbine.csv not found, default created. ' ...
                'Units: [s, ms^-1]']}];
        end
    case 'Constant'
        try
            Wind.Vel = readmatrix('WindVelConstant.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'WindVelConstant.csv',1,[],8,[])
            loadDataWarnings = [loadDataWarnings{:},...
                {['WindVelConstant.csv not found, default created. ' ...
                'Unit: [ms^-1]']}];
        end
    case 'ZOH_wErrorCov'
        Wind.Vel.Init = readmatrix('WindVelConstant.csv');
        VelCov = readmatrix('WindVelCovariance.csv');
        nT = size(turbProp.Pos,1);
        [~,Wind.Vel.CholSig] = readCovMatrix(VelCov,nT,'WindVel');
    case 'RW_with_Mean'
        Wind.Vel.Init = readmatrix('WindVelConstant.csv');
        VelCov = readmatrix('WindVelCovariance.csv');
        nT = size(turbProp.Pos,1);
        [~,Wind.Vel.CholSig] = readCovMatrix(VelCov,nT,'WindVel');
        % MeanPull \in [0,1] scales the influence of the difference between
        % the mean and the current state. A smaller value leads to bigger
        % value "detours" from the mean.
        Wind.Vel.MeanPull = 1;
    case 'Interpolation_wErrorCov'
        Wind.Vel.Data = readmatrix('WindVel.csv');
        VelCov = readmatrix('WindVelCovariance.csv');
        nT = size(turbProp.Pos,1);
        [~,Wind.Vel.CholSig] = readCovMatrix(VelCov,nT,'WindVel');
    case 'InterpTurbine_wErrorCov'
        Wind.Vel.Data = readmatrix('WindVelTurbine.csv');
        VelCov = readmatrix('WindVelCovariance.csv');
        nT = size(turbProp.Pos,1);
        [~,Wind.Vel.CholSig] = readCovMatrix(VelCov,nT,'WindVel');
    case 'Constant_wErrorCov'
        Wind.Vel.Data = readmatrix('WindVelConstant.csv');
        VelCov = readmatrix('WindVelCovariance.csv');
        nT = size(turbProp.Pos,1);
        [~,Wind.Vel.CholSig] = readCovMatrix(VelCov,nT,'WindVel');
    otherwise
        error(['Method for wind velocity ' Wind.Input_Vel ' unknown.'])
end
% ============= Direction ============= 
switch Wind.Input_Dir
    case 'Interpolation'
        try
            Wind.Dir = readmatrix('WindDir.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'WindDir.csv',2,size(turbProp.Pos,1),...
                [Sim.StartTime,270],[Sim.EndTime,225])
            loadDataWarnings = [loadDataWarnings{:},...
                {['WindDir.csv not found, default created. ' ...
                'Units: [s, deg]']}];
        end
        % =============================================================== %
        % FIX transition from e.g. 10 to 350 deg and the other way around
        % =============================================================== %
    case 'InterpTurbine'
        try
            Wind.Dir = readmatrix('WindDirTurbine.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'WindDirTurbine.csv',3,size(turbProp.Pos,1),...
                [Sim.StartTime,270],[Sim.EndTime,225])
            loadDataWarnings = [loadDataWarnings{:},...
                {['WindDirTurbine.csv not found, default created. ' ...
                'Units: [s, deg]']}];
        end
        % =============================================================== %
        % FIX transition from e.g. 10 to 350 deg and the other way around
        % =============================================================== %
    case 'Constant'
        try
            Wind.Dir = readmatrix('WindDirConstant.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'WindDirConstant.csv',1,[],270,[])
            loadDataWarnings = [loadDataWarnings{:},...
                {['WindDirConstant.csv not found, default created. ' ...
                'Unit: [deg]']}];
        end
    case 'Interpolation_wErrorCov'
        Wind.Dir.Data = readmatrix('WindDir.csv');
        DirCov = readmatrix('WindDirCovariance.csv');
        nT = size(turbProp.Pos,1);
        [~,Wind.Dir.CholSig] = readCovMatrix(DirCov,nT,'WindDir');
        % =============================================================== %
        % FIX transition from e.g. 10 to 350 deg and the other way around
        % =============================================================== %
    case 'InterpTurbine_wErrorCov'
        Wind.Dir.Data = readmatrix('WindDirTurbine.csv');
        DirCov = readmatrix('WindDirCovariance.csv');
        nT = size(turbProp.Pos,1);
        [~,Wind.Dir.CholSig] = readCovMatrix(DirCov,nT,'WindDir');
        % =============================================================== %
        % FIX transition from e.g. 10 to 350 deg and the other way around
        % =============================================================== %
    case 'Constant_wErrorCov'
        Wind.Dir.Data = readmatrix('WindDirConstant.csv');
        DirCov = readmatrix('WindDirCovariance.csv');
        nT = size(turbProp.Pos,1);
        [~,Wind.Dir.CholSig] = readCovMatrix(DirCov,nT,'WindDir');
        % =============================================================== %
        % FIX transition from e.g. 10 to 350 deg and the other way around
        % =============================================================== %
    case 'RW_with_Mean'
        Wind.Dir.Init = readmatrix('WindDirConstant.csv');
        DirCov = readmatrix('WindDirCovariance.csv');
        nT = size(turbProp.Pos,1);
        [~,Wind.Dir.CholSig] = readCovMatrix(DirCov,nT,'WindDir');
        % MeanPull \in [0,1] scales the influence of the difference between
        % the mean and the current state. A smaller value leads to bigger
        % value "detours" from the mean.
        Wind.Dir.MeanPull = 1;
        % =============================================================== %
        % FIX transition from e.g. 10 to 350 deg and the other way around
        % =============================================================== %
    otherwise
        error(['Method for wind direction ' Wind.Input_Dir ' unknown.'])
end
% ============= TI ============= 
switch Wind.Input_TI
    case 'Interpolation'
        try
            Wind.TI = readmatrix('WindTI.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'WindTI.csv',2,[],...
                [Sim.StartTime,0.06],[Sim.EndTime,0.10])
            loadDataWarnings = [loadDataWarnings{:},...
                {['WindTI.csv not found, default created. ' ...
                'Units: [s, percent]']}];
        end
    case 'InterpTurbine'
        Wind.TI = readmatrix('WindTITurbine.csv');
        try
            Wind.TI = readmatrix('WindTITurbine.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'WindTITurbine.csv',3,size(turbProp.Pos,1),...
                [Sim.StartTime,0.06],[Sim.EndTime,0.10])
            loadDataWarnings = [loadDataWarnings{:},...
                {['WindTITurbine.csv not found, default created. ' ...
                'Units: [s, percent]']}];
        end
    case 'Constant'
        try
            Wind.TI = readmatrix('WindTIConstant.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'WindTIConstant.csv',1,[],0.06,[])
            loadDataWarnings = [loadDataWarnings{:},...
                {['WindTIConstant.csv not found, default created. ' ...
                'Unit: [percent]']}];
        end
    otherwise
        error(['Method for turbulence intensity ' Wind.Input_TI ' unknown.'])
end
% ============= Shear ============= 
switch Wind.Input_Shr
    case 'PowerLaw'
        alpha = readmatrix('WindShearPowerLaw.csv');
        Wind.Shear.alpha = alpha;
        Wind.Shear.z0 = 1;
        clear alpha
    case 'Interpolation'
        Wind.Shear = readmatrix('WindShearProfile.csv');
    case 'LogLaw'
        Wind.Shear = readmatrix('WindShearLogLaw.csv');
    otherwise
        error(['Method for wind shear ' Wind.Input_Shr ' unknown.'])
end

%% Prepare turbine struct
% Turbine positions
T.posBase   = turbProp.Pos;
T.nT        = size(turbProp.Pos,1);
t_data      = getTurbineData(turbProp.Type);
T.posNac    = t_data.NacPos;
T.D         = t_data.D;

% States / OP chains
states      = States();

if isequal(paramFLORIDyn.twf_model, 'heterogeneous')
    % OP orientation is used to derive TWF, needs to be added.
    states.WF_names{end+1} = 'OP_ori';
    states.WF              = length(states.WF_names);
elseif ~isequal(paramFLORIDyn.twf_model, 'homogeneous')
    error(['The TWF method "' paramFLORIDyn.twf_model '" is not known.' ...
        ' Use "homogeneous" or "heterogeneous" instead.'])
end

%   OP coordinates
T.States_OP = zeros(paramFLORIDyn.n_op * T.nT, states.OP);
T.Names_OP  = states.OP_names;
%   Turbine states
T.States_T  = zeros(paramFLORIDyn.n_op * T.nT, states.Turbine);
T.Names_T   = states.T_names;
%   Wind field states
T.States_WF = zeros(paramFLORIDyn.n_op * T.nT, states.WF);
T.Names_WF  = states.WF_names;
%   Starting ID of all OP chains and number of OPs
T.StartI    = 1:paramFLORIDyn.n_op:paramFLORIDyn.n_op * T.nT;
T.nOP       = paramFLORIDyn.n_op;
%   Interaction matrix
%       1 -> no influence
%       0 -> Reduces wind speed to 0.
T.red_arr   = ones(T.nT,T.nT);

% Check if upstream dependency is assigned and assign if it isn't
if ~isfield(paramFLORIDyn,'deltaUW')
    paramFLORIDyn.deltaUW = paramFLORIDyn.deltaDW;
end
%% Load Control parameters
switch Con.Yaw%     = 'Interploation'; % 'Interploation','SOWFA'
    case 'Constant'
        try
            Con.YawData = readmatrix('Control_YawConstant.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'Control_YawConstant.csv',1,[],270,[])
            loadDataWarnings = [loadDataWarnings{:},...
                {['Control_YawConstant.csv not found, default created. ' ...
                'Unit: [deg]']}];
        end
    case 'InterpTurbine'
        try
            Con.YawData = readmatrix('Control_YawInterpolation.csv');
        catch
            generateDemoCSV([Sim.PathToSim 'Data' filesep],...
                'Control_YawInterpolation.csv',3,size(turbProp.Pos,1),...
                [Sim.StartTime,270],[Sim.EndTime,250])
            loadDataWarnings = [loadDataWarnings{:},...
                {['Control_YawInterpolation.csv not found, default created. ' ...
                'Unit: [deg]']}];
        end
    case 'SOWFA'
        nacelleYaw  = importSOWFAFile([Sim.PathToSim 'Data' filesep...
            'SOWFA_nacelleYaw.csv']);
        Con.YawData  = condenseSOWFAYaw(...
            [nacelleYaw(1:T.nT:end,2),...
            reshape(nacelleYaw(:,3),T.nT,[])']);
    otherwise
        error(['Yawing method ' Con.Yaw ' unknown.'])
end

if ~isfield(Con,'tanhYaw'); Con.tanhYaw = false; end

%% Initialize OP states
[T.States_OP, T.States_T, T.States_WF] = InitStates(...
    T,Wind,turbProp.Init_States,paramFLORIS,Sim);

%% Sim
Sim.nSimSteps = length(Sim.StartTime:Sim.TimeStep:Sim.EndTime);
paramFLORIS.RotorPoints = Sim.RotorPoints;

%% Vis
if Vis.FlowField.Plot.Online
    Vis.Film.MovFileEffU = [Sim.PathToSim filesep 'Results' filesep ...
        'EffectiveWindSpeed.avi'];
    Vis.Film.FrmFileEffU(Sim.nSimSteps) = struct('cdata',[],'colormap',[]);
    Vis.Film.InProgress = true;
else
    Vis.Film.InProgress = false;
end

if Vis.FlowField.Error.Online
    dirSnap = dir(Vis.FlowField.Error.ValidationPath);
    Vis.FlowField.Error.Steps = zeros(length(dirSnap)-2,1);
    for i = 3:length(dirSnap)
        try
            Vis.FlowField.Error.Steps(i-2) = str2num(dirSnap(i).name);
        catch
            Vis.FlowField.Error.Steps(i-2) = -1;
        end
    end
    %Vis.FlowField.Error.Steps = str2num(dirSap
end

%% Check if there were issues loading data
if ~isempty(loadDataWarnings)
    for iW = 1:length(loadDataWarnings)
        warning(loadDataWarnings{iW});
    end
    error(['Data could not be read as intended, default files have ' ... 
        'been generated (see warnings) - please change the data in the' ...
        ' generated files to what is intended.'])
end
end

function generateDemoCSV(path,name,type,nT,startV,endV)
switch type
    case 1
        % Constant value
        %   start: single value
        writematrix(startV,[path name])
    case 2
        % Interpolation
        %   start: [time, value]
        %   end:   [time, value]
        writematrix([startV;endV],[path name])
    case 3
        % Turbine individual interpolation
        %   start: [time, value]
        %   end:   [time, value]
        A = [startV(2);endV(2)];
        T = [startV(1);endV(1)];
        writematrix([T,repmat(A,1,nT)],[path name])
    otherwise
        error(['prepareSimulation -> generateDemoCSV: Type ' ...
            num2str(type) ' not defined.'])
end
    
end