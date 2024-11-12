pathEnKF = [Sim.PathToSim 'Results' filesep 'FlowFieldRaw'];
pathOrig = [Sim.PathToSim 'Results' filesep 'FlowFieldRaw_Original'];

napshots = dir([Sim.PathToSim 'Results' filesep 'FlowFieldRaw' filesep 'Dir']);
meas = {'Dir','Vel','Eff'};

%% 
fieldRes  = Vis.FlowField.Res;
fieldLims = Vis.FlowField.Lims;
xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);

[X,Y] = meshgrid(xAx,yAx);

%%

f = figure;

clf;
gap_HW    = [.06 .04];
marg_h_LU = [.06 .06];
marg_w_LR = [.02 .02];
ha = tight_subplot(3,3,gap_HW,marg_h_LU,marg_w_LR);
f = gcf;

for step = 1:length(napshots)-2
    timestep = napshots(step+2).name;
    OPs = readmatrix([pathOrig filesep 'OP' filesep timestep]);
    OPsE = readmatrix([pathEnKF filesep 'OP' filesep timestep]);
    for i = 1:3
        M_O = readmatrix([pathOrig filesep meas{i} filesep timestep]);
        M_E = readmatrix([pathEnKF filesep meas{i} filesep timestep]);

        % Calc limits
        limits_O = [floor(min(M_O,[],"all")) ceil(max(M_O,[],"all"))];
        limits_E = [floor(min(M_E,[],"all")) ceil(max(M_E,[],"all"))];
        limits_O = [min(limits_O(1),limits_E(1)) max(limits_O(2),limits_E(2))];
        maxErr = ceil(max(abs(M_O-M_E),[],"all"));
        limits_R = [-maxErr maxErr];

        %% FLORIDyn
        axes(ha(1+(i-1)*3))
        contourf(X,Y,M_O,10,'EdgeColor','none');
        hold on
        scatter(OPs(:,1),OPs(:,2),5,0.9*ones(length(OPs),3),'filled')
        scatter(OPs(1:10:end,1),OPs(1:10:end,2),15,...
            0.9*ones(length(OPs(1:10:end,1)),3),'filled')
        hold off
        axis equal

        c = colorbar(ha(1+(i-1)*3),"westoutside");

        colormap(ha(1+(i-1)*3),viridis(11))
        caxis(limits_O)

        xlim([fieldLims(1,1),fieldLims(2,1)])
        ylim([fieldLims(1,2),fieldLims(2,2)])
        title(['Original state at t = ' timestep(1:end-4) ' s']);
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        %% EnKF
        axes(ha(2+(i-1)*3))
        contourf(X,Y,M_E,10,'EdgeColor','none');
        hold on
        scatter(OPsE(:,1),OPsE(:,2),5,zeros(length(OPs),3),'filled')
        scatter(OPsE(1:10:end,1),OPsE(1:10:end,2),15,...
            zeros(length(OPs(1:10:end,1)),3),'filled')
        hold off
        colormap(ha(2+(i-1)*3),viridis(11))
        c = colorbar(ha(2+(i-1)*3),"westoutside");
        caxis(limits_O)
        axis equal
        xlim([fieldLims(1,1),fieldLims(2,1)])
        ylim([fieldLims(1,2),fieldLims(2,2)])
        title(['EnKF state at t = ' timestep(1:end-4) ' s']);
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        %% Error
        axes(ha(3+(i-1)*3))
        contourf(X,Y,M_O-M_E,10,'EdgeColor','none');
        hold on
        scatter(OPsE(:,1),OPsE(:,2),5,zeros(length(OPs),3),'filled')
        scatter(OPsE(1:10:end,1),OPsE(1:10:end,2),15,...
            zeros(length(OPs(1:10:end,1)),3),'filled')

        scatter(OPs(:,1),OPs(:,2),5,0.9*ones(length(OPs),3),'filled')
        scatter(OPs(1:10:end,1),OPs(1:10:end,2),15,...
            0.9*ones(length(OPs(1:10:end,1)),3),'filled')
        hold off
        colormap(ha(3+(i-1)*3),flipud(RdBu(11)))
        c = colorbar;
        title(['Error at t = ' timestep(1:end-4) ' s']);
        caxis(limits_R)
        axis equal
        xlim([fieldLims(1,1),fieldLims(2,1)])
        ylim([fieldLims(1,2),fieldLims(2,2)])
        set(gca,'XTick',[])
        set(gca,'YTick',[])
    end
    if not(isfolder([pathEnKF filesep 'Comparison']))
        mkdir([pathEnKF filesep 'Comparison']); 
    end
    saveas(f,[pathEnKF filesep 'Comparison' filesep timestep(1:end-4) '.png'])
end

