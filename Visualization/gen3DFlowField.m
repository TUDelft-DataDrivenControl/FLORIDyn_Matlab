function [Vis] = gen3DFlowField(T,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS)
%PLOTFLOWFIELD Generate full flow field plot
% exportgraphics(gcf,'9T_flow.pdf','ContentType','vector')
%% Preallocate field
nM = countMeasurements(Vis.FlowField.Data);
nM = 3;
fieldLims = Vis.FlowField.Lims;
fieldRes  = Vis.FlowField3D.Res;
xAx = fieldLims(1,1):fieldRes(1):fieldLims(2,1);
yAx = fieldLims(1,2):fieldRes(2):fieldLims(2,2);
zAx = fieldLims(1,3):fieldRes(3):fieldLims(2,3);

[X,Y,H] = meshgrid(xAx,yAx,zAx);

clear fieldRes 

%% Get data

if Vis.FlowField3D.Parallel
    [Zred,ZambI,Zeff] = getMeasurementsPar2(X,Y,H,nM,T,paramFLORIS,Wind,Vis);
else
    [Zred,ZambI,Zeff] = getMeasurements(X,Y,H,nM,T,paramFLORIS,Wind,Vis);
end

[~,~,~] = mkdir([Sim.PathToSim 'Results']);
%% Generate & save plots
% figure
% contourf(X,Y,Z(:,:,1),40,'LineColor','none')
% axis equal
% colorbar
% title('Wind Speed Reduction')

%% Added turbulence
if Vis.FlowField.Data.AddedTurb
    Ztmp = sqrt(ZambI);
    Ztmp(Ztmp>1) = 0;
    genVTKField(X,Y,H,Ztmp,...
        [Sim.PathToSim 'Results\addedTurbulence' ],'AddedTurbulence');
end
%% Effective Wind Speed
if Vis.FlowField.Data.WindSpeedEff
    genVTKField(X,Y,H,Zeff,...
        [Sim.PathToSim 'Results\effWindSpeed' ],'EffectiveWindSpeed');
end
end
%exportgraphics(gcf,'9T_HorFlow_t300.pdf','ContentType','vector')
%%
function [Zred,ZambI,Zeff] = getMeasurements(X,Y,H,nM,T,paramFLORIS,Wind,Vis)
% Wrapper function to disguise the grid points as turbines with one rotor
% point and experience almost the same calculations as the rotor points in
% the simulation.
%   Single thread version
%sizeX = size(X);
Zred  = zeros(size(X));
Zeff  = zeros(size(X));
ZambI = zeros(size(X));

%Z = zeros(sizeX(1),sizeX(2),sizeX(3),nM);
%parfor iGP = 1:length(X(:))
for iGP = 1:length(X(:))
    xGP = X(iGP);
    yGP = Y(iGP);
    zGP = H(iGP);
    
    GPdep = {1:T.nT};


    % Create T equivalent to get interpolated OPs
    GP.dep          = cell(T.nT,1);
    GP.dep(1)       = GPdep(1);
    GP.StartI       = T.StartI;
    GP.nT           = 1;
    GP.States_OP    = T.States_OP;
    GP.posBase      = [xGP,yGP,0];
    GP.nOP          = T.nOP;
    GP.intOPs = interpolateOPs(GP);
    
    GP.posNac = [0,0,zGP];
    GP.States_WF = T.States_WF;
    GP.States_T = T.States_T;
    GP.D = [T.D;0];
    
    [tmpM,~] = setUpTmpWFAndRun(GP,paramFLORIS,Wind);
    
    Zred(iGP)  = tmpM(1);
    ZambI(iGP) = tmpM(2);
    Zeff(iGP)  = tmpM(3);
%     [rw,cl,cl2] = ind2sub(size(X),iGP);
%     Z(rw,cl,cl2,1:3) = tmpM;
end
end

function [Zred,ZambI,Zeff] = getMeasurementsPar2(X,Y,H,nM,T,paramFLORIS,Wind,Vis)
% Wrapper function to disguise the grid points as turbines with one rotor
% point and experience almost the same calculations as the rotor points in
% the simulation.
%   Multi thread version, requires parallel computing toolbox
%   For one plot, the start-up time of the parallel toolbox can be longer
%   than the saved computational time.
Zred  = zeros(size(X));
Zeff  = zeros(size(X));
ZambI = zeros(size(X));

%parfor iGP = 1:length(X(:))
parfor iGP = 1:length(X(:))
    xGP = X(iGP);
    yGP = Y(iGP);
    zGP = H(iGP);
    
    GPdep = {1:T.nT};

    GP = struct();
    % Create T equivalent to get interpolated OPs
    GP.dep          = cell(T.nT,1);
    GP.dep(1)       = GPdep(1);
    GP.StartI       = T.StartI;
    GP.nT           = 1;
    GP.States_OP    = T.States_OP;
    GP.posBase      = [xGP,yGP,0];
    GP.nOP          = T.nOP;
    GP.intOPs = interpolateOPs(GP);
    
    GP.posNac = [0,0,zGP];
    GP.States_WF = T.States_WF;
    GP.States_T = T.States_T;
    GP.D = [T.D;0];
    
    [tmpM,~] = setUpTmpWFAndRun(GP,paramFLORIS,Wind);
    
    Zred(iGP)  = tmpM(1);
    ZambI(iGP) = tmpM(2);
    Zeff(iGP)  = tmpM(3);
end

%Z = cat(3,Zred,ZambI,Zeff);
end

function Z = getMeasurementsDelaunay(T)
DT = delaunay(T.States_OP(:,1),T.States_OP(:,2));
trisurf(DT,T.States_OP(:,1),T.States_OP(:,2),T.States_WF(:,1),...
    'EdgeColor','none','facecolor', 'interp')
end
%%
function Z = Wrap_Z(Z,tmpM,sizeX)
[rw,cl] = ind2sub(sizeX,iGP);
    tmpM1 = sub2ind(sizeZ,rw,cl,1);
    tmpM2 = sub2ind(sizeZ,rw,cl,2);
    tmpM3 = sub2ind(sizeZ,rw,cl,3);
    
    Z(tmpM1) = tmpM(1);
end
%%
function n = countMeasurements(VisFlowFieldData)
n = 0;
if VisFlowFieldData.WindSpeedFree
    n = n+1;
end
if VisFlowFieldData.WindSpeedEff
    n = n+1;
end
if VisFlowFieldData.Reduction
    n = n+1;
end
if VisFlowFieldData.AddedTurb
    n = n+1;
end
if VisFlowFieldData.WindDirection
    n = n+1;
end
end
% Vis.FlowField.Plot.Post             = false;
% Vis.FlowField.Plot.Online           = false;
% Vis.FlowField.Plot.OnlineFreq       = 1; % every n-th time step
% Vis.FlowField.Data.WindSpeedFree    = false;
% Vis.FlowField.Data.WindSpeedEff     = true;
% Vis.FlowField.Data.Reduction        = false;
% Vis.FlowField.Data.AddedTurb        = false;
% Vis.FlowField.Data.WindDirection    = false;
% Vis.FlowField.Lims = [0 0 0; 2100 2100 400]; % [min x,y,z; max x,y,z] in m
% Vis.FlowField.Res = 10; % in m
% Vis.FlowField.Save = false;