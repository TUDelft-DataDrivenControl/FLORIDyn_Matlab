function errorVTK(path2VTK,filename,T,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS)
%ERRORVTK means to create a VTK file with the errors of FLORIDyn on
%comparison to a validation vtk file. Will use the etracted cellCenters to
%find the matching values in FLORIDyn.

% CURRENTLY ONLY FOR 2D HORIZONTAL SNAPSHOTS

%% Extract data
if strcmp(path2VTK(end-4:end),'.vtk')
    [~,cellCenters,cellData] = importVTK(path2VTK);
else
    [~,cellCenters,cellData] = importVTK([path2VTK '.vtk']);
end

Validation_UmeanAbs = sqrt(sum(cellData.^2,2));
H = mean(cellCenters(:,3));
%% Get FLORIDyn data
X = cellCenters(:,1);
Y = cellCenters(:,2);
if Vis.FlowField.Plot.parallel
    F_U = getMeasurementsPar2(X,Y,3,H,T,paramFLORIS,Wind,Vis);
else
    F_U = getMeasurements(X,Y,3,H,T,paramFLORIS,Wind,Vis);
end
%% Calculate error
errU = Validation_UmeanAbs - F_U;
%% Save
%p=[reshape(X,prod(size(X)),1),reshape(Y,prod(size(Y)),1)];
p=[X,Y];
t=delaunay(p);

if not(isfolder([Sim.PathToSim '\Results\ErrorVTK']))
    mkdir([Sim.PathToSim '\Results\ErrorVTK'])
end

writeVTK([Sim.PathToSim '\Results\ErrorVTK\' filename],t,p,errU)

%% Plot
persistent extE muSig

if isempty(extE)
    extE = 0;
end
extE = max(abs([errU(:);extE]));

fplot = figure(99);
subplot(1,2,1)
trisurf(t,X,Y,errU,'EdgeColor','none')
xlabel('West-East in m')
ylabel('South-North in m')
%colormap(TUDelft_Div_rwb(100))
colormap(RdBu)
c = colorbar;
caxis([-8,8])
c.Label.String = 'Val.-FLORIDyn in ms^{-1}';
axis equal
view(2)
title(['t = ' filename(1:end-9) ' s'])

subplot(1,2,2)
if isempty(muSig)
    muSig = [mean(errU),std(errU)];
else
    muSig = [muSig;mean(errU),std(errU)];
    d_mu = gradient(muSig(:,1)) ;
    d_sig = gradient(muSig(:,2)) ;
end



%plot(muSig(:,1),muSig(:,2),'LineWidth',0.5,'Color',[235,114,70]./255)
if size(muSig,1)>1
    quiver(muSig(1:end-1,1),muSig(1:end-1,2),d_mu(1:end-1),d_sig(1:end-1),...
        0,'Color',[235,114,70]./255)
    hold on
    scatter(muSig(:,1),muSig(:,2),8,'filled','MarkerEdgeColor','none',...
              'MarkerFaceColor',[0,113,136]./255)
    hold off
else
    scatter(muSig(:,1),muSig(:,2),8,'filled','MarkerEdgeColor','none',...
              'MarkerFaceColor',[0,113,136]./255)
end
ylim([0,2])
xlim([-2,2])
grid on
title('Mean and std of the error')
ylabel('\sigma in ms^{-1}')
xlabel('\mu in ms^{-1}')

% if firstPlot
%     plot([0,0],[-1,2],'--k','LineWidth',1)
%     xlabel('Error [ms^{-1}]')
%     ylabel('Probability')
% end
% hold on
% y       = linspace(-10,10);
% mu      = mean(errU);
% sigma   = std(errU);
% f       = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',0.5,'Color',[235,114,70]./255)
% ylim([0,1])
% xlim([-10,10])
% grid on
% title('Error PDF')
fplot.Position(3:4) = [1147, 420];
drawnow
saveas(fplot,[Sim.PathToSim '\Results\ErrorVTK\Error' filename(1:end-9) '.png'])
end

function U = getMeasurements(X,Y,~,zh,T,paramFLORIS,Wind,Vis)
% Wrapper function to disguise the grid points as turbines with one rotor
% point and experience almost the same calculations as the rotor points in
% the simulation.
%   Single thread version
nPoints = length(X);
U = zeros(nPoints,1);
%parfor iGP = 1:length(X(:))
for iGP = 1:nPoints
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
    U(iGP) = tmpM(3);
%     [rw,cl] = ind2sub(size(X),iGP);
%     U(rw,cl,1:3) = tmpM;
end
end

function U = getMeasurementsPar2(X,Y,~,zh,T,paramFLORIS,Wind,~)
% Wrapper function to disguise the grid points as turbines with one rotor
% point and experience almost the same calculations as the rotor points in
% the simulation.
%   Multi thread version, requires parallel computing toolbox
%   For one plot, the start-up time of the parallel toolbox can be longer
%   than the saved computational time.
%Zred  = zeros(size(X));
Zeff  = zeros(size(X));
%ZambI = zeros(size(X));

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
    
    %Zred(iGP)  = tmpM(1);
    %ZambI(iGP) = tmpM(2);
    Zeff(iGP)  = tmpM(3);
end
U = Zeff;
%Z = cat(3,Zred,ZambI,Zeff);
end