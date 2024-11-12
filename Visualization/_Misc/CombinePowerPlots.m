


%%
if exist([Sim.PathToSim 'Data\SOWFA_generatorPower.csv'],'file')
    powSOWFA = importSOWFAFile([Sim.PathToSim 'Data\SOWFA_generatorPower.csv']);
    powSOWFA(:,3) = powSOWFA(:,3)./(10^6*paramFLORIS.airDen);
    pSOW = true;
elseif exist([Sim.PathToSim 'Data/FC_generatorPower.csv'],'file')
    powSOWFAtmp = readmatrix('FC_generatorPower.csv');
    
    pow_tmp = powSOWFAtmp(:,2:end)';
    t_tmp = repmat(powSOWFAtmp(:,1)',T.nT,1);
    
    powSOWFA = [zeros(size(t_tmp(:))),t_tmp(:),pow_tmp(:)];
    
    pSOW = true;
else
    pSOW = false;
end

if Vis.Msmnts.Data.Probes
    addpath([Sim.PathToSim '/Data'])
    Vis.Msmnts.Data.Probes = exist('U.csv','file')>0;
    if ~Vis.Msmnts.Data.Probes
        warning('Probe file U.csv not found, no probe data plotted.')
    end
end
%% Subtract staring time 
if Vis.SubtractOffset
    timeFDyn = M.("Time [s]") - M.("Time [s]")(1)- Vis.Msmnts.Data.TimeOffset;
    if pSOW
        timeSOW = powSOWFA(:,2) - powSOWFA(1,2)- Vis.Msmnts.Data.TimeOffset;
    end
    
else
    timeFDyn = M.("Time [s]");
    if pSOW
        timeSOW = powSOWFA(:,2);
    end
end
%%

lay = getLayout(T.nT);
            
if pSOW
    y_lim = [max(min([M.("Power generated [MW]");powSOWFA(:,3)]),0),...
        max([M.("Power generated [MW]");powSOWFA(:,3)])];
else
    y_lim = [max(min(M.("Power generated [MW]")),0),...
        max(M.("Power generated [MW]"))];
end
y_lim = y_lim + [-0.1 0.1]*max((y_lim(2)-y_lim(1)),0.5);
y_lim(1) = max(y_lim(1),0);
f = figure;

c = inferno(3);
for iT = 1:T.nT
    subplot(lay(1),lay(2),iT)
    px = [0 1 1 0]*300+600;
    py = [0 0 1 1]*12;
    patch(px,py,[1,1,1]*0.8,'FaceAlpha',0.3,'EdgeColor','none')
    hold on
    plot(...
        timeFDyn(iT:T.nT:end),...
        M.("Power generated [MW]")(iT:T.nT:end),...
        'LineWidth',1.5,'Color',c(1,:));
    if pSOW
        %hold on
        plot(...
            timeFDyn(iT:T.nT:end),...
            MInI.("Power generated [MW]")(iT:T.nT:end),...
            'LineWidth',1.5,'Color',c(2,:));
        plot(...
            timeSOW(iT:T.nT:end),...
            powSOWFA(iT:T.nT:end,3),...
            'LineWidth',1.5,'Color',c(3,:));
        hold off
        legend('','FLORIDyn w/o I&I','FLORIDyn w I&I','Validation','location','best')
    end
    grid on
    title(['Turbine ' num2str(iT-1)])
    xlim([max(timeFDyn(1),0),timeFDyn(end)])
    ylim(y_lim)
    xlabel('Time [s]')
    ylabel('Power generated [MW]')
end

%exportgraphics(gcf,'name.pdf','ContentType','vector')