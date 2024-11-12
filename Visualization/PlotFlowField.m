function [Vis] = PlotFlowField(T,~,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS,SimTime)
%PLOTFLOWFIELD Generate full flow field plot
% exportgraphics(gcf,'9T_flow.pdf','ContentType','vector')
%% Preallocate field
nM = countMeasurements(Vis.FlowField.Data);
nM = 3;
fieldLims = Vis.FlowField.Lims;
fieldRes  = Vis.FlowField.Res;
xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);

[X,Y] = meshgrid(xAx,yAx);

clear fieldRes 

%% Get data

if Vis.FlowField.Plot.parallel
    Z = getMeasurementsPar2(X,Y,nM,T.posNac(1,3),T,paramFLORIS,Wind,Vis);
else
    Z = getMeasurements(X,Y,nM,T.posNac(1,3),T,paramFLORIS,Wind,Vis);
end

%% Generate & save plots
% figure
% contourf(X,Y,Z(:,:,1),40,'LineColor','none')
% axis equal
% colorbar
% title('Wind Speed Reduction')

%% Added turbulence
if Vis.FlowField.Data.AddedTurb
    figure
    Ztmp = sqrt(Z(:,:,2));
    Ztmp(Ztmp>1) = 0;
    hold on
    cLeveles = 10;
    contourf(X,Y,Ztmp,cLeveles,'LineColor','none')
    % Plot Rotors
    for i_T = 1:length(T.D)
        % Get start and end of the turbine rotor
        yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2) - ...
            T.States_T(T.StartI(i_T),2));
        %yaw = deg2rad(T.States_T(T.StartI(i_T),2));
        rot_pos = ...
            [cos(yaw), -sin(yaw);...
            sin(yaw), cos(yaw)] * ...
            [0,0;T.D(i_T)/2,-T.D(i_T)/2];
        rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
        plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'k','LineWidth',3);
    end
    hold off
    axis equal
    c = colorbar;
    c.Label.String ='Added turbulence';
    colormap(viridis(cLeveles))
    xlabel('West-East (m)')
    ylabel('South-North (m)')
    title('Added turbulence')
    f.Units               = 'centimeters';
    f.Position(3)         = 16.1; % line width
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.04))
end
%% Effective Wind Speed
if Vis.FlowField.Data.WindSpeedEff
    if Vis.FlowField.Plot.Online
        f = figure(101);
        clf
    else
        f = figure;
    end
    set(gcf,'color','w');
    contourf(X,Y,Z(:,:,3),40,'LineColor','none')
    hold on
    % Plot Rotors
    for i_T = 1:length(T.D)
        % Get start and end of the turbine rotor
        %yaw = deg2rad(T.States_T(T.StartI(i_T),2));
        yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
            T.States_T(T.StartI(i_T),2));
        rot_pos = ...
            [cos(yaw), -sin(yaw);...
            sin(yaw), cos(yaw)] * ...
            [0,0;T.D(i_T)/2,-T.D(i_T)/2];
        rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
        plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'k','LineWidth',3);
    end
    if Vis.FlowField.Plot.OPs
        scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
            'filled','MarkerFaceColor',[1 1 1])
        scatter(T.States_OP(1:10:end,1),T.States_OP(1:10:end,2),15,...
            'filled','MarkerFaceColor',[1 1 1])
    end
    axis equal
    xlim([fieldLims(1,1),fieldLims(2,1)])
    ylim([fieldLims(1,2),fieldLims(2,2)])
    
    hold off
    
    c = colorbar;
    c.Label.String ='Wind speed (ms^{-1})';
    colormap(inferno(1000))
    xlabel('West-East (m)')
    ylabel('South-North (m)')
    if Vis.Film.InProgress
        title(['Effective Wind Speed, t = ' num2str(SimTime) ' s'])
    else
        title('Effective Wind Speed')
    end
    f.Units               = 'centimeters';
    f.Position(3)         = 16.1; % line width
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.04))
    c.Limits = [2, 10];                                     % HARDCODED FOR PLOTTING!!!
    clim([2, 10]);                                          % HARDCODED FOR PLOTTING!!!
    drawnow
    
    if and(Vis.FlowField.Plot.Online,Vis.Film.InProgress)
        c.Limits = Vis.Film.LimFileEffU;
        set(gca, 'Clim', Vis.Film.LimFileEffU)
        caxis(Vis.Film.LimFileEffU)
        colormap(viridis(1000))
        %colormap(TUDelft_Seq_wyr(1000));
        % Save fram for movie
        Vis.Film.FrmFileEffU(Sim.SimStep) = getframe(f);
        if Sim.SimStep == Sim.nSimSteps
            Vis.Film.InProgress = false;
        end
    end
end
end
%exportgraphics(gcf,'9T_HorFlow_t300.pdf','ContentType','vector')
%%
function Z = getMeasurements(X,Y,nM,zh,T,paramFLORIS,Wind,Vis)
% Wrapper function to disguise the grid points as turbines with one rotor
% point and experience almost the same calculations as the rotor points in
% the simulation.
%   Single thread version
sizeX = size(X);
Z = zeros(sizeX(1),sizeX(2),nM);
%parfor iGP = 1:length(X(:))
for iGP = 1:length(X(:))
    xGP = X(iGP);
    yGP = Y(iGP);
    
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
    
    GP.posNac = [0,0,zh];
    GP.States_WF = T.States_WF;
    GP.States_T = T.States_T;
    GP.D = [T.D;0];
    
    [tmpM,~] = setUpTmpWFAndRun(GP,paramFLORIS,Wind);
    
    [rw,cl] = ind2sub(size(X),iGP);
    Z(rw,cl,1:3) = tmpM;
end
end

function Z = getMeasurementsPar2(X,Y,nM,zh,T,paramFLORIS,Wind,Vis)
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
    
    GP.posNac = [0,0,zh];
    GP.States_WF = T.States_WF;
    GP.States_T = T.States_T;
    GP.D = [T.D;0];
    
    [tmpM,~] = setUpTmpWFAndRun(GP,paramFLORIS,Wind);
    
    Zred(iGP)  = tmpM(1);
    ZambI(iGP) = tmpM(2);
    Zeff(iGP)  = tmpM(3);
end

Z = cat(3,Zred,ZambI,Zeff);
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
%exportgraphics(gcf,'name.pdf','ContentType','vector')
%c = colorbar;
%c.Limits = [0,14];
%set(gca, 'Clim', [0, 14])
%
%