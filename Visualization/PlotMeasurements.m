% Copyright (C) <2024>, M Becker
%
% List of the contributors to the development of FLORIDyn: see LICENSE file.
% Description and complete License: see LICENSE file.
	
% This program (FLORIDyn) is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program (see COPYING file).  If not, see <https://www.gnu.org/licenses/>.
% ======================================================================= %
% Updated: 16. Dez. 2024, M. Becker
% ======================================================================= %

function [] = PlotMeasurements(T,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS)
%Plot all requested measurements

%% Load SOWFA if existent
if exist([Sim.PathToSim 'Data' filesep 'SOWFA_generatorPower.csv'],'file')
    powSOWFA = importSOWFAFile([Sim.PathToSim 'Data' filesep 'SOWFA_generatorPower.csv']);
    powSOWFA(:,3) = powSOWFA(:,3)./(10^6*paramFLORIS.airDen);
    pSOW = true;
elseif exist([Sim.PathToSim 'Data' filesep 'FC_generatorPower.csv'],'file')
    powSOWFAtmp = readmatrix('FC_generatorPower.csv');
    
    pow_tmp = powSOWFAtmp(:,2:end)';
    t_tmp = repmat(powSOWFAtmp(:,1)',T.nT,1);
    
    powSOWFA = [zeros(size(t_tmp(:))),t_tmp(:),pow_tmp(:)];
    
    pSOW = true;
else
    pSOW = false;
end

if Vis.Msmnts.Data.Probes
    addpath([Sim.PathToSim filesep 'Data'])
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
%% Power generated
if Vis.Msmnts.Data.Power
    
    switch Vis.Msmnts.Layout
        case 'separated'
            lay = getLayout(T.nT);
            
            if pSOW
                y_lim = [min([M.("Power generated [MW]");powSOWFA(:,3)]),...
                    max([M.("Power generated [MW]");powSOWFA(:,3)])];
            else
                y_lim = [min(M.("Power generated [MW]")),...
                    max(M.("Power generated [MW]"))];
            end
            y_lim = y_lim + [-0.1 0.1]*max((y_lim(2)-y_lim(1)),0.5);
            
            f = figure;
            if pSOW
                c = inferno(2);
            else
                c = inferno(1);
            end
            for iT = 1:T.nT
                subplot(lay(1),lay(2),iT)
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Power generated [MW]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c(1,:));
                if pSOW
                    hold on
                    plot(...
                        timeSOW(iT:T.nT:end),...
                        powSOWFA(iT:T.nT:end,3),...
                        'LineWidth',1,'Color',c(2,:));
                    hold off
                    legend('FLORIDyn','Validation','location','best')
                end
                grid on
                title(['Turbine ' num2str(iT-1)])
                xlim([max(timeFDyn(1),0),timeFDyn(end)])
                ylim(y_lim)
                xlabel('Time [s]')
                ylabel('Power generated [MW]')
            end
            %TUcolors;
        otherwise
            f = figure;
            c = inferno(T.nT);
            hold on
            for iT = 1:T.nT
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Power generated [MW]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c(iT,:));
                if pSOW
                    plot(...
                        timeSOW(iT:T.nT:end),...
                        powSOWFA(iT:T.nT:end,3),...
                        'LineWidth',1,'Color',c(iT,:),'LineStyle','-.');
                end
            end
            grid on
            title('Power generated [MW]')
            xlim([max(timeFDyn(1),0),timeFDyn(end)])
            xlabel('Time [s]')
            ylabel('Power generated [MW]')
            hold off
    end
end
%% Power Difference between SOWFA and FLORIDyn
if and(Vis.Msmnts.Data.PowerDiffSOWFA,pSOW)
    switch Vis.Msmnts.Layout
        case 'separated'
            lay = getLayout(T.nT);
            
%             y_lim = [min(M.("Power generated [MW]")),...
%                 max(M.("Power generated [MW]"))];
%             y_lim = y_lim + [-0.1 0.1]*max((y_lim(2)-y_lim(1)),0.5);
            
            f = figure;
            c = inferno(1);
            for iT = 1:T.nT
                deltaP = interp1(timeSOW(iT:T.nT:end),...
                    powSOWFA(iT:T.nT:end,3),...
                    timeFDyn(iT:T.nT:end))...
                    - M.("Power generated [MW]")(iT:T.nT:end);
                subplot(lay(1),lay(2),iT)
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    deltaP,...
                    'LineWidth',2,'Color',c);
                
                grid on
                title(['Turbine ' num2str(iT-1)])
                xlim([max(timeFDyn(1),0),timeFDyn(end)])
                %                     ylim(y_lim)
                xlabel('Time [s]')
                ylabel('Delta Power [MW]')
            end
        otherwise
            f = figure;
            c = inferno(T.nT);
            hold on
            for iT = 1:T.nT
                deltaP = interp1(timeSOW(iT:T.nT:end),...
                    powSOWFA(iT:T.nT:end,3),...
                    timeFDyn(iT:T.nT:end))...
                    - M.("Power generated [MW]")(iT:T.nT:end);
%                 plot(...
%                     timeFDyn(iT:T.nT:end),...
%                     deltaP,...
%                     'LineWidth',2,'Color',c(iT,:));
                scatter(...
                    M.("Effective Wind Speed [ms^-1]")(iT:T.nT:end),...
                    deltaP,10,c(iT,:),'filled')%...
                %'LineWidth',2,'Color',c(iT,:));
            end
            grid on
            title('Power difference SOWFA - FLORIDyn [MW]')
            xlim([max(timeFDyn(1),0),timeFDyn(end)])
            xlabel('Eff. Wind Speed [s]')
            ylabel('Delta Power [MW]')
            hold off
    end
else
    if(Vis.Msmnts.Data.PowerDiffSOWFA)
        warning(['Can not find power generated in the Validation simulation,'...
            ' will not plot the difference.'])
    end
end

%% Effective Wind Speed
if Vis.Msmnts.Data.EffWindSpeed
    switch Vis.Msmnts.Layout
        case 'separated'
            lay = getLayout(T.nT);
            
            y_lim = [min(M.("Effective Wind Speed [ms^-1]")),...
                max(M.("Effective Wind Speed [ms^-1]"))];
            y_lim = y_lim + [-0.1 0.1]*max((y_lim(2)-y_lim(1)),0.5);
            
            f = figure;
            c = inferno(1);
            for iT = 1:T.nT
                subplot(lay(1),lay(2),iT)
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Effective Wind Speed [ms^-1]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c);
                grid on
                title(['Turbine ' num2str(iT-1)])
                xlim([max(timeFDyn(1),0),timeFDyn(end)])
                ylim(y_lim)
                xlabel('Time [s]')
                ylabel('Eff. wind speed [ms^{-1}]')
            end
        otherwise
            f = figure;
            c = inferno(T.nT);
            hold on
            for iT = 1:T.nT
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Effective Wind Speed [ms^-1]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c(iT,:));
            end
            grid on
            title('Effective Wind Speed [ms^{-1}]')
            xlim([max(timeFDyn(1),0),timeFDyn(end)])
            xlabel('Time [s]')
            ylabel('Eff. wind speed [ms^{-1}]')
            hold off
    end
end

%% Added turbulence
if Vis.Msmnts.Data.AddedTurb
    switch Vis.Msmnts.Layout
        case 'separated'
            lay = getLayout(T.nT);
            
            y_lim = [min(M.("Added turbulence [%]")),...
                max(M.("Added turbulence [%]"))];
            y_lim = y_lim + [-0.1 0.1]*max((y_lim(2)-y_lim(1)),0.5);
            
            f = figure;
            c = inferno(1);
            for iT = 1:T.nT
                subplot(lay(1),lay(2),iT)
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Added turbulence [%]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c);
                grid on
                title(['Turbine ' num2str(iT-1)])
                xlim([max(timeFDyn(1),0),timeFDyn(end)])
                ylim(y_lim)
                xlabel('Time [s]')
                ylabel('Added turbulence [%]')
            end
        otherwise
            f = figure;
            c = inferno(T.nT);
            hold on
            for iT = 1:T.nT
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Added turbulence [%]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c(iT,:));
            end
            grid on
            title('Added turbulence [%]')
            xlim([max(timeFDyn(1),0),timeFDyn(end)])
            xlabel('Time [s]')
            ylabel('Added turbulence [%]')
            hold off
    end
end

%% Foreign Reduction
if Vis.Msmnts.Data.ForeignRed
    switch Vis.Msmnts.Layout
        case 'separated'
            lay = getLayout(T.nT);
            
            y_lim = [min(M.("Foreign Reduction [%]")),...
                max(M.("Foreign Reduction [%]"))];
            y_lim = y_lim + [-0.1 0.1]*max((y_lim(2)-y_lim(1)),0.5);
            
            f = figure;
            c = inferno(1);
            for iT = 1:T.nT
                subplot(lay(1),lay(2),iT)
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Foreign Reduction [%]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c);
                grid on
                title(['Turbine ' num2str(iT-1)])
                xlim([max(timeFDyn(1),0),timeFDyn(end)])
                ylim(y_lim)
                xlabel('Time [s]')
                ylabel('Foreign Reduction [%]')
            end
        otherwise
            f = figure;
            c = inferno(T.nT);
            hold on
            for iT = 1:T.nT
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Foreign Reduction [%]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c(iT,:));
            end
            grid on
            title('Foreign Reduction [%]')
            xlim([max(timeFDyn(1),0),timeFDyn(end)])
            xlabel('Time [s]')
            ylabel('Foreign Reduction [%]')
            hold off
    end
end

%% Free wind speed
if Vis.Msmnts.Data.FreeWindSpeed
    switch Vis.Msmnts.Layout
        case 'separated'
            lay = getLayout(T.nT);
            
            y_lim = [min(M.("Free Wind Speed [ms^-1]")),...
                max(M.("Free Wind Speed [ms^-1]"))];
            y_lim = y_lim + [-0.1 0.1]*max((y_lim(2)-y_lim(1)),0.5);
            
            f = figure;
            c = inferno(1);
            for iT = 1:T.nT
                subplot(lay(1),lay(2),iT)
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Free Wind Speed [ms^-1]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c);
                grid on
                title(['Turbine ' num2str(iT-1)])
                xlim([max(timeFDyn(1),0),timeFDyn(end)])
                ylim(y_lim)
                xlabel('Time [s]')
                ylabel('Free wind speed [ms^{-1}]')
            end
        otherwise
            f = figure;
            c = inferno(T.nT);
            hold on
            for iT = 1:T.nT
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Free Wind Speed [ms^-1]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c(iT,:));
            end
            grid on
            title('Free Wind Speed [ms^{-1}]')
            xlim([max(timeFDyn(1),0),timeFDyn(end)])
            xlabel('Time [s]')
            ylabel('Free wind speed [ms^{-1}]')
            hold off
    end
end

%% Free and effective wind speed
if Vis.Msmnts.Data.EffFreeWindSpeed
    switch Vis.Msmnts.Layout
        case 'separated'
            lay = getLayout(T.nT);
            
            y_lim = [min(...
                [M.("Free Wind Speed [ms^-1]");...
                 M.("Effective Wind Speed [ms^-1]")]),...
                max(...
                [M.("Free Wind Speed [ms^-1]");...
                 M.("Effective Wind Speed [ms^-1]")])];
            y_lim = y_lim + [-0.1 0.1]*max((y_lim(2)-y_lim(1)),0.5);
            
            f = figure;
            if Vis.Msmnts.Data.Probes
                c = inferno(3);
                U = readmatrix('U.csv');
                Udiff = zeros(size(M.("Free Wind Speed [ms^-1]"),1),2);
                Udiff(:,1) = timeFDyn;
            else
                c = inferno(2);
            end
            
            for iT = 1:T.nT
                subplot(lay(1),lay(2),iT)
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Free Wind Speed [ms^-1]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c(1,:));
                hold on
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Effective Wind Speed [ms^-1]")(iT:T.nT:end),...
                    'LineWidth',1,'Color',c(2,:));
                
                if Vis.Msmnts.Data.Probes
                    plot(...
                        U(:,1)-U(1,1),...
                        U(:,1+iT),...
                        'LineWidth',2,'Color',c(3,:));
                end
                
                hold off
                grid on
                title(['Turbine ' num2str(iT-1)])
                xlim([max(timeFDyn(1),0),timeFDyn(end)])
                ylim(y_lim)
                xlabel('Time [s]')
                ylabel('Wind speed [ms^{-1}]')
                if Vis.Msmnts.Data.Probes
                    legend('Free (I&I \div FLORIDyn)','Effective (I&I)',...
                        'True free','Location','northwest')
                    Udiff(iT:T.nT:end,2) = ...
                        interp1(U(:,1)-U(1,1),U(:,1+iT),...
                        timeFDyn(iT:T.nT:end)) - ...
                        M.("Free Wind Speed [ms^-1]")(iT:T.nT:end);
                else
                    legend('Free (I&I \div FLORIDyn)','Effective (I&I)',...
                        'Location','northwest')
                end
                
            end
            
            % Plotting probe data
            if Vis.Msmnts.Data.Probes
                figure
                
                y_lim = [min(Udiff(:,2)),...
                    max(Udiff(:,2))];
                y_lim = y_lim + [-0.1 0.1]*max((y_lim(2)-y_lim(1)),0.5);
                
                for iT = 1:T.nT
                    subplot(lay(1),lay(2),iT)
                    plot(...
                        timeFDyn(iT:T.nT:end),...
                        Udiff(iT:T.nT:end,2),...
                        'LineWidth',2,'Color',c(1,:));
                    grid on
                    title(['Turbine ' num2str(iT-1)])
                    xlim([max(timeFDyn(1),0),timeFDyn(end)])
                    ylim(y_lim)
                    xlabel('Time [s]')
                    ylabel('delta Wind speed [ms^{-1}]')
                    legend('diff true free and est. free',...
                        'Location','northwest')
                end
            end
        otherwise
            f = figure;
            c = inferno(T.nT);
            hold on
            for iT = 1:T.nT
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Free Wind Speed [ms^-1]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c(iT,:));
                plot(...
                    timeFDyn(iT:T.nT:end),...
                    M.("Effective Wind Speed [ms^-1]")(iT:T.nT:end),...
                    'LineWidth',2,'Color',c(iT,:),'LineStyle','--');
            end
            grid on
            title('Free and effective wind speed [ms^{-1}]')
            xlim([max(timeFDyn(1),0),timeFDyn(end)])
            xlabel('Time [s]')
            ylabel('Wind speed [ms^{-1}]')
            hold off
    end
end
end
%exportgraphics(gcf,'name.pdf','ContentType','vector')
