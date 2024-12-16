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


classdef EnKF_Visualizer
    %ENKF_VISUALIZER Function to facilitate the plotting of EnKF related
    %data
    
    properties
        Settings
        Figs
        AxesHandles
        tmpLims
    end
    
    methods
        function obj = EnKF_Visualizer(Settings)
            %ENKF_VISUALIZER Construct an instance of this class
            %   Detailed explanation goes here
            obj.Settings = Settings;
            
            % Figs holds links to all figure objects
            %   0 + _ -> state & correction during simulation (live)
            %   1 -> velocity 
            %   2 -> direction
            %   3 -> Turbunlence Intensity
            %
            %   3 + _ -> Localization (once)
            %   1 -> Velocity 
            %   2 -> Direction
            %   3 -> Turbunlence Intensity
            %
            %   7 -> Undefined
            %
            %   7 + _ -> State and variance
            %   1 -> Velocity 
            %   2 -> Direction
            %   3 -> Turbunlence Intensity
            %
            %   11 -> Quiver OPs
            %
            %   12 -> Flow field
            %
            %   12 + _ -> Correction of each state by each turbine
            %   1 -> Velocity 
            %   2 -> Direction
            %   3 -> Turbunlence Intensity
            %   
            %   15 + _ -> combined measurement
            %   1 -> Foreign Reduction
            %   2 -> Added turbulence
            %   3 -> Effective Wind Speed
            %   4 -> Free Wind Speed
            %   5 -> Power Generated
            %   6 -> Wind Direction
            %
            %   22 -> Error plot
            %   
            %   22 + _ -> Quiver correction plot
            %   1 -> Velocity
            %   2 -> Wind Dir
            %   3 -> Added turbulence
            %
            %   25 + _ -> flow field
            %   1 -> Effective wind speed
            %   2 -> Direction
            %   3 -> Added turbulence
            %
            %   29 + _ -> Weighted State and variance
            %   1 -> Velocity 
            %   2 -> Direction
            %   3 -> Turbunlence Intensity
            %
            %   32 -> Weighted wind speed field and variance w SOWFA

            obj.Figs = cell(32,1);
            obj.AxesHandles = cell(32,1);
            
            obj.tmpLims = [0,1];
        end
        
        function obj = plotK_combined(obj, nT, EnKF, Vis, K, diff, iE, iState)
            %% plotK_combined creates a live figure with all OPs
            % from the first ensembe and in which direction they are
            % corrected. The left subplot shows the correction, the right
            % the state value.
            % iState
            %   1 -> Velocity
            %   2 -> Direction
            %   3 -> Turbulence intensity
            if iE>1
                return
            end
            %plotK visualizes where the measurement of each turbine has an
            %influence on the states due to K
            if isempty(obj.Figs{iState});obj.Figs{iState} = figure; end
            try
                figure(obj.Figs{iState})
            catch
                obj.Figs{iState} = figure;
            end
            clf;
            fieldLims = Vis.FlowField.Lims;
            
            %% Get Data
            x = EnKF.States_OP(:,1); % mean(EnKF.States_OP(:,1:EnKF.nStatesOP:end),2);
            y = EnKF.States_OP(:,2); %mean(EnKF.States_OP(:,2:EnKF.nStatesOP:end),2);
            c = K*diff;
            max_c = max(abs(c));
            
            %% Plotting
            ax(1) = subplot(1,2,1);
            scatter(x,y,[],c,'filled')
            colormap(ax(1),RdBu)
            caxis([-max_c,max_c])
            grid on
            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])
            colb = colorbar;
            colb.Label.String = 'Correction';
            
            
            ax(2) = subplot(1,2,2);
            switch iState
                case 1
                    scatter(x,y,[],EnKF.States.Vel(:,iE),'filled')
                case 2
                    scatter(x,y,[],EnKF.States.Dir(:,iE),'filled')
                case 3
                    scatter(x,y,[],EnKF.States.TI(:,iE),'filled')
            end
            grid on
            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])
            colormap(ax(2),inferno);
            colb = colorbar;
            colb.Label.String = 'State Values';
            
            obj.Figs{iState}.Position(1:4) = [200, 200, 1380, 600];
            pause(0.1)
        end
        
        function obj = plotK_Localization(obj, T, Vis, coLength, iState)
            %% Creates a Contour plot for all turbines to show how far their
            % influence reaches, called once.
            if isempty(obj.Figs{3+iState});obj.Figs{3+iState} = figure; end
            try
                figure(obj.Figs{3+iState})
            catch
                obj.Figs{3+iState} = figure;
            end
            clf;
            
            % Create subplot
            l = getLayout(T.nT);
            fieldLims = Vis.FlowField.Lims;
            fieldRes  = Vis.FlowField.Res;
            xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
            yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);
            
            [X,Y] = meshgrid(xAx,yAx);
            for iT = 1:T.nT
                subplot(l(1),l(2),iT)
                distGP = sqrt((T.posBase(iT,1) - X(:)).^2 +...
                    (T.posBase(iT,2) - Y(:)).^2);
                LocCov = GaspariAndCohn1999(coLength,distGP);
                contourf(X,Y,reshape(LocCov,size(X)))
                colormap(viridis)
                axis equal
                xlim([fieldLims(1,1),fieldLims(2,1)])
                ylim([fieldLims(1,2),fieldLims(2,2)])
            end
            
        end
        
        function obj = plotK_PowerGenerated(obj, EnKF,T)
            %% Plots the power generated by all ensembles, split by turbine
            %   Only recommended to get an idea of the raw data
            if isempty(obj.Figs{7});obj.Figs{7} = figure; end
            try
                figure(obj.Figs{7})
            catch
                obj.Figs{7} = figure;
            end
            clf;
            
            l = getLayout(T.nT);
            maxP = 0;
            minP = 100;
            for iE = 1:EnKF.nE
                for iT = 1:T.nT
                    VisTime = EnKF.M{iE}.("Time [s]");
                    VisPow  = EnKF.M{iE}.("Power generated [MW]");
                    maxP = max(maxP,max(VisPow));
                    minP = min(minP,min(VisPow));
                    subplot(l(1),l(2),iT)
                    hold on
                    plot(VisTime(iT:T.nT:end)-VisTime(1),...
                        VisPow(iT:T.nT:end),...
                        'Color',EnKF.Colors(mod(iE,11)+1,:))
                    hold off
                    grid on
                end
            end
            
            for iT = 1:T.nT
                subplot(l(1),l(2),iT)
                offS = (maxP-minP)*0.1;
                ylim([minP-offS,maxP+offS])
                ylabel('Power generated (MW)')
                xlabel('Time (s)')
                title(['Turbine' num2str(iT-1)]) 
            end
        end
        
        function obj = plotK_CombinedMeasurement(obj, EnKF, T, Sim, paramFLORIS, iM,indT)
            %% Plots the mean power generated with +-3,2,1 standard
            % deviations
            if isempty(obj.Figs{15+iM});obj.Figs{15+iM} = figure; end
            try
                figure(obj.Figs{15+iM})
            catch
                obj.Figs{15+iM} = figure;
            end
            l = getLayout(length(indT));
            gap_HW    = [.05 .05];
            marg_h_LU = [.06 .01];
            marg_w_LR = [.04 .02];
            ha = tight_subplot(l(1),l(2),gap_HW,marg_h_LU,marg_w_LR);
            f = gcf;
            f.Position = [520 190 830 650];
            
            % Get combined measurements and extract mean power + std
            M_comb = cat(1,EnKF.M{:});
            VisTime = EnKF.M{1}.('Time [s]');
            VisTime = VisTime(1:T.nT:end)-VisTime(1);
            
            powSOWFA = false;
            switch iM
                case 1
                    Data_raw = reshape(...
                        M_comb.("Foreign Reduction [%]"),[],EnKF.nE);
                    yLabel = 'Foreign Reduction (%)';
                case 2
                    Data_raw = reshape(...
                        M_comb.("Added turbulence [%]"),[],EnKF.nE);
                    yLabel = 'Added turbulence (%)';
                case 3
                    Data_raw = reshape(...
                        M_comb.("Effective Wind Speed [ms^-1]"),[],EnKF.nE);
                    yLabel = 'Effective Wind Speed (ms^-1)';
                case 4
                    Data_raw = reshape(...
                        M_comb.("Free Wind Speed [ms^-1]"),[],EnKF.nE);
                    yLabel = 'Free Wind Speed (ms^-1)';
                case 5
                    Data_raw = reshape(...
                        M_comb.("Power generated [MW]"),[],EnKF.nE);
                    yLabel = 'Power generated (MW)';
                    
                    if exist([Sim.PathToSim 'Data' filesep 'SOWFA_generatorPower.csv'],'file')
                        powSOWFA = importSOWFAFile([Sim.PathToSim 'Data' filesep 'SOWFA_generatorPower.csv']);
                        if isempty(powSOWFA)
                            powSOWFA = csvread([Sim.PathToSim 'Data' filesep 'SOWFA_generatorPower.csv']);
                        end
                        powSOWFA(:,3) = powSOWFA(:,3)./(10^6*paramFLORIS.airDen);
                        timeSOW = powSOWFA(:,2) - powSOWFA(1,2);
                        pSOW = true;
                    elseif exist([Sim.PathToSim 'Data' filesep 'FC_generatorPower.csv'],'file')
                        powSOWFAtmp = readmatrix('FC_generatorPower.csv');
                        
                        pow_tmp = powSOWFAtmp(:,2:end)';
                        t_tmp = repmat(powSOWFAtmp(:,1)',T.nT,1);
                        
                        powSOWFA = [zeros(size(t_tmp(:))),t_tmp(:),pow_tmp(:)];
                        timeSOW = powSOWFA(:,2) - powSOWFA(1,2);
                        pSOW = true;
                    elseif exist([Sim.PathToSim 'Data' filesep 'FLORIDyn_generatorPower.csv'],'file')
                        powSOWFA = readmatrix('FLORIDyn_generatorPower.csv');
                        timeSOW = powSOWFA(:,2) - powSOWFA(1,2);
                        pSOW = true;
                    else
                        pSOW = false;
                    end
                case 6
                    error('Wind direction not supported yet as measurement.')
                otherwise
                    error(['Measurement ' num2str(iM) ' not known'])
            end
            
            Data_Mean = mean(Data_raw,2);
            Data_Stdv = std(Data_raw,[],2);
            
            maxD = max(Data_Mean + 3*Data_Stdv);
            minD = min(Data_Mean - 3*Data_Stdv);
            if((maxD-minD) == 0)
                offS = 1;
            else
                offS = (maxD-minD)*0.1;
            end
            
            % Plot for all turbines individually
            for i = 1:length(indT)
                iT = indT(i);
                axes(ha(i))
                if iM == 5
                    if pSOW
                        plot(...
                            timeSOW(iT:T.nT:end),...
                            powSOWFA(iT:T.nT:end,3),'-',...
                            'LineWidth',1,'Color',EnKF.Colors(1,:));
                        hold on
                    end
                end
                fill([VisTime' fliplr(VisTime')],...
                    [(Data_Mean(iT:T.nT:end) + 3*Data_Stdv(iT:T.nT:end))',...
                    fliplr((Data_Mean(iT:T.nT:end) - 3*Data_Stdv(iT:T.nT:end))')],...
                    EnKF.Colors(5,:),'FaceAlpha',0.1,'EdgeColor','none')%,...
                    %'EdgeColor',sqrt(EnKF.Colors(5,:)),'LineWidth',0.1)
                hold on
                fill([VisTime' fliplr(VisTime')],...
                    [(Data_Mean(iT:T.nT:end) + 2*Data_Stdv(iT:T.nT:end))',...
                    fliplr((Data_Mean(iT:T.nT:end) - 2*Data_Stdv(iT:T.nT:end))')],...
                    EnKF.Colors(5,:),'FaceAlpha',0.2,'EdgeColor','none')%,...
                    %'EdgeColor',sqrt(EnKF.Colors(5,:)),'LineWidth',0.1)
                fill([VisTime' fliplr(VisTime')],...
                    [(Data_Mean(iT:T.nT:end) + 1*Data_Stdv(iT:T.nT:end))',...
                    fliplr((Data_Mean(iT:T.nT:end) - 1*Data_Stdv(iT:T.nT:end))')],...
                    EnKF.Colors(5,:),'FaceAlpha',0.3,'EdgeColor','none')%,...
                    %'EdgeColor',sqrt(EnKF.Colors(5,:)),'LineWidth',0.1)
                plot(VisTime,...
                    Data_Mean(iT:T.nT:end),...
                    'Color',EnKF.Colors(5,:),'LineWidth',1.2)
                
                hold off
                grid on
                
                ylim([max(minD-offS,0),maxD+offS])
                xlim([min(VisTime),max(VisTime)])
                if i == floor(l(1)/2)*l(2)+1
                    ylabel(yLabel)
                elseif i == (l(1)-1)*l(2)+round(l(2)/2)
                    xlabel('Time (s)')
                end
            end
            
        end
        
        function obj = plotK_StateAndVariance(obj, EnKF, T, iState, val, Vis)
            %% Plots the average states with +- 3,2,1 sigma bounds 
            % and the absolute variance in one plot. Can be a bit cluttered
            if isempty(obj.Figs{7+iState});obj.Figs{7+iState} = figure; end
            try
                figure(obj.Figs{7+iState})
            catch
                obj.Figs{7+iState} = figure;
            end
            clf;
            
            switch iState
                case 1
                    Sts = 1:size(EnKF.States.Vel,1);
                    Std = std(EnKF.States.Vel,[],2);
                    Mean = mean(EnKF.States.Vel,2);
                    yLabel = 'Wind speed (m/s)';
                    yyLabel = 'Std (m/s)';
                    yLimsS = [6.5,9.5];
                    yLimsV = [0,0.8];
                case 2
                    Sts = 1:size(EnKF.States.Dir,1);
                    Std = std(EnKF.States.Dir,[],2);
                    Mean = mean(EnKF.States.Dir,2);
                    yLabel = 'Wind direction (deg)';
                    yyLabel = 'Std (deg)';
                    yLimsS = [val-5,val+5];
                    yLimsV = [0,3];
                case 3
                    Sts = 1:size(EnKF.States.TI,1);
                    Std = std(EnKF.States.TI,[],2);
                    Mean = mean(EnKF.States.TI,2);
                    yLabel = 'Amb. Turb. Intensity (%)';
                    yyLabel = 'Std (%)';
                otherwise
                    error('State index not known, use a value between 1 & 3')
            end
            l = getLayout(T.nT);
            gap_HW    = [.04 .07];
            marg_h_LU = [.06 .06];
            marg_w_LR = [.06 .06];
            ha = tight_subplot(l(1),l(2),gap_HW,marg_h_LU,marg_w_LR);
            f = gcf;
            f.Position = [520 190 830 650];

            trueOPpos = ...
                reshape( ...
                mean(...
                reshape(...
                EnKF.States_OP,...
                length(T.Names_OP)*T.nT*T.nOP,[])...
                ,2),...
                T.nT*T.nOP,[]);

%             yLimsS = [min(Mean-3*Std) max(Mean+3*Std)]; 
%             yLimsS = yLimsS + [-.1 .1]*(yLimsS(2)-yLimsS(1));
%             yLimsV = [min(Std) max(Std)]; 
%             yLimsV = yLimsV + [-.1 .1]*(yLimsV(2)-yLimsV(1));
            for iT = 1:T.nT
                axes(ha(iT))
                colororder([EnKF.Colors(5,:);EnKF.Colors(1,:)])
                si = (iT-1)*T.nOP + 1;
                ei = iT*T.nOP;
                
                % 
                withIn = and(trueOPpos(si:ei,1) < Vis.FlowField.Lims(2,1),...
                    trueOPpos(si:ei,2) < Vis.FlowField.Lims(2,2));
                indBorder = find(withIn==0,1);
                
                yyaxis left
                try
                    fill([Sts(si+indBorder-1),Sts(ei),...
                        Sts(ei), Sts(si+indBorder-1)],...
                        [yLimsS(1),yLimsS(1),yLimsS(2),yLimsS(2)],...
                        EnKF.Colors(1,:),'FaceAlpha',0.29,'EdgeColor','none')
                    hold on
                    fill([Sts(si:ei), fliplr(Sts(si:ei))],...
                        [(Mean(si:ei) + 3*Std(si:ei));...
                        flipud((Mean(si:ei) - 3*Std(si:ei)))],...
                        EnKF.Colors(5,:),'FaceAlpha',0.10,'EdgeColor','none')
                catch
                    fill([Sts(si:ei), fliplr(Sts(si:ei))],...
                        [(Mean(si:ei) + 3*Std(si:ei));...
                        flipud((Mean(si:ei) - 3*Std(si:ei)))],...
                        EnKF.Colors(5,:),'FaceAlpha',0.10,'EdgeColor','none')
                    hold on
                end
                fill([Sts(si:ei), fliplr(Sts(si:ei))],...
                    [(Mean(si:ei) + 2*Std(si:ei));...
                    flipud((Mean(si:ei) - 2*Std(si:ei)))],...
                    EnKF.Colors(5,:),'FaceAlpha',0.15,'EdgeColor','none')
                fill([Sts(si:ei), fliplr(Sts(si:ei))],...
                    [(Mean(si:ei) + 1*Std(si:ei));...
                    flipud((Mean(si:ei) - 1*Std(si:ei)))],...
                    EnKF.Colors(5,:),'FaceAlpha',0.20,'EdgeColor','none')
                plot(Sts(si:ei),Mean(si:ei),'-',...
                    'Color',EnKF.Colors(5,:),'LineWidth',1.2)
                
                if nargin >= 5
                    plot([Sts(si), Sts(ei)],[val, val],':',...
                        'Color',0.8*EnKF.Colors(5,:),'LineWidth',1.2)
                end
                hold off
                ax = gca;
                ax.XTick = si:10:ei;
                grid on
                
                ylim(yLimsS)
                %ylim([2,12])
                if iT == floor(l(1)/2)*l(2)+1
                    ylabel(yLabel)
                end
                
                yyaxis right
                plot(Sts(si:ei),Std(si:ei),'-.',...
                    'Color',EnKF.Colors(1,:),'LineWidth',0.4)
                
                % Add moving average
                hold on
                fw = 10; B = 1/fw*ones(fw,1);
%                 inOut = ones(fw-1,1)*Std(si);
%                 initCond = filtic(B,1,inOut,inOut);
%                 plot(Sts(si:ei),filter(B,1,Std(si:ei),initCond),'-',...
%                     'Color',EnKF.Colors(1,:),'LineWidth',1.2)
                plot(Sts(si:ei),filtfilt(B,1,Std(si:ei)),'-',...
                    'Color',EnKF.Colors(1,:),'LineWidth',1.2)
                hold off
                ylim(yLimsV)
                %ylim([0,3])
                title(['Turbine ' num2str(iT-1)])
                
                xlim([(iT-1)*T.nOP+1 iT*T.nOP])
                if iT == floor(l(1)/2)*l(2)+l(2)
                    ylabel(yyLabel)
                elseif iT == (l(1)-1)*l(2)+round(l(2)/2)
                    xlabel('State')
                end
            end
        end
        
        function obj = plotK_WeightedStateAndVariance(obj, EnKF, T, Vis, Sim, iState,val)
            %% Plots the average states with +- 3,2,1 sigma bounds 
            % and the absolute variance in one plot. Can be a bit cluttered
            if isempty(obj.Figs{7+iState});obj.Figs{7+iState} = figure; end
            try
                set(0,'CurrentFigure',figure(obj.Figs{7+iState}));
            catch
                obj.Figs{7+iState} = figure;
            end
            clf;

            [trueOPpos, Mean, ~, Std] = ...
                EnKF_WeightedVariance(EnKF, T, Vis, Sim, iState);
            dist2closestTurbine = min(sqrt(...
                (T.posBase(:,1)'-trueOPpos(:,1)).^2 + ...
                (T.posBase(:,2)'-trueOPpos(:,2)).^2),[],2);

            switch iState
                case 1
                    Sts = 1:size(EnKF.States.Vel,1);
%                     Std = std(EnKF.States.Vel,[],2);
%                     Mean = mean(EnKF.States.Vel,2);
                    yLabel = 'Wind speed (m/s)';
                    yyLabel = 'Std (m/s)';
                case 2
                    Sts = 1:size(EnKF.States.Dir,1);
%                     Std = std(EnKF.States.Dir,[],2);
%                     Mean = mean(EnKF.States.Dir,2);
                    yLabel = 'Wind direction (deg)';
                    yyLabel = 'Std (deg)';
                case 3
                    Sts = 1:size(EnKF.States.TI,1);
%                     Std = std(EnKF.States.TI,[],2);
%                     Mean = mean(EnKF.States.TI,2);
                    yLabel = 'Amb. Turb. Intensity (%)';
                    yyLabel = 'Std (%)';
                otherwise
                    error('State index not known, use a value between 1 & 3')
            end
            l = getLayout(T.nT);
            gap_HW    = [.04 .07];
            marg_h_LU = [.06 .06];
            marg_w_LR = [.06 .06];
            ha = tight_subplot(l(1),l(2),gap_HW,marg_h_LU,marg_w_LR);
            f = gcf;
            f.Position = [520 190 830 650];
            
            
            yLimsS = [min(Mean-3*Std) max(Mean+3*Std)]; 
            yLimsS = yLimsS + [-.1 .1]*(yLimsS(2)-yLimsS(1));
            yLimsV = [min(Std) max(Std)]; 
            yLimsV = yLimsV + [-.1 .1]*(yLimsV(2)-yLimsV(1));
            for iT = 1:T.nT
                axes(ha(iT))
                colororder([EnKF.Colors(5,:);EnKF.Colors(1,:)])
                si = (iT-1)*T.nOP + 1;
                ei = iT*T.nOP;
                yyaxis left
                fill([Sts(si:ei), fliplr(Sts(si:ei))],...
                    [(Mean(si:ei) + 3*Std(si:ei));...
                    flipud((Mean(si:ei) - 3*Std(si:ei)))],...
                    EnKF.Colors(5,:),'FaceAlpha',0.10,'EdgeColor','none')
                hold on
                fill([Sts(si:ei), fliplr(Sts(si:ei))],...
                    [(Mean(si:ei) + 2*Std(si:ei));...
                    flipud((Mean(si:ei) - 2*Std(si:ei)))],...
                    EnKF.Colors(5,:),'FaceAlpha',0.15,'EdgeColor','none')
                fill([Sts(si:ei), fliplr(Sts(si:ei))],...
                    [(Mean(si:ei) + 1*Std(si:ei));...
                    flipud((Mean(si:ei) - 1*Std(si:ei)))],...
                    EnKF.Colors(5,:),'FaceAlpha',0.20,'EdgeColor','none')
                plot(Sts(si:ei),Mean(si:ei),'-',...
                    'Color',EnKF.Colors(5,:),'LineWidth',1.2)
                
                if nargin == 7
                    plot([Sts(si), Sts(ei)],[val, val],':',...
                        'Color',0.8*EnKF.Colors(5,:),'LineWidth',1.2)
                end
                hold off
                grid on
                
                ylim(yLimsS)
                %ylim([2,12])
                if iT == floor(l(1)/2)*l(2)+1
                    ylabel(yLabel)
                end
                
                yyaxis right
                plot(Sts(si:ei),Std(si:ei),'-.',...
                    'Color',EnKF.Colors(1,:),'LineWidth',0.4)
                
                % Add moving average
                hold on
                fw = 10; B = 1/fw*ones(fw,1);
%                 inOut = ones(fw-1,1)*Std(si);
%                 initCond = filtic(B,1,inOut,inOut);
%                 plot(Sts(si:ei),filter(B,1,Std(si:ei),initCond),'-',...
%                     'Color',EnKF.Colors(1,:),'LineWidth',1.2)
                plot(Sts(si:ei),filtfilt(B,1,Std(si:ei)),'-',...
                    'Color',EnKF.Colors(1,:),'LineWidth',1.2)
                hold off
                ylim(yLimsV)
                %ylim([0,3])
                title(['Turbine ' num2str(iT-1)])
                
                xlim([(iT-1)*T.nOP+1 iT*T.nOP])
                if iT == floor(l(1)/2)*l(2)+l(2)
                    ylabel(yyLabel)
                elseif iT == (l(1)-1)*l(2)+round(l(2)/2)
                    xlabel('State')
                end
            end
        end


        function obj = plotK_ArchOPs(obj, EnKF, T, Vis)
            %% Creates a flow field plot with an arch for every op which
            % indicates wind speed and direction. Not useful for large flow
            % fields as the wind speed is plotted with 1m:1m/s scale
            if isempty(obj.Figs{11});obj.Figs{11} = figure; end
            try
                figure(obj.Figs{11})
            catch
                obj.Figs{11} = figure;
            end
            clf;
            
            % See which values are varied
            varU_Phi = false(2,1);
            if EnKF.Vel.Correct
                varU_Phi(1) = true;
                
                dataU = EnKF.States.Vel_mean;
                stdU = std(EnKF.States.Vel,[],2);
            else
                dataU = T.States_WF(:,1);
            end
            
            if EnKF.Dir.Correct
                varU_Phi(2) = true;
                dataPhi = EnKF.States.Dir_mean;
                stdPhi = std(EnKF.States.Dir,[],2);
            else
                dataPhi = T.States_WF(:,2);
            end
            OP_meanX = EnKF.States_OP(:,1:EnKF.nStatesOP:end);
            OP_meanY = EnKF.States_OP(:,2:EnKF.nStatesOP:end);
            %
            if and(varU_Phi(1),varU_Phi(2))
                % Plot variable wind speed and direction
                hold on
                for iOP = 1:length(dataU)
                    phiW = angSOWFA2world(dataPhi(iOP));
                    [x,y] = getArch(...
                        dataU(iOP),stdU(iOP),...
                        phiW,deg2rad(stdPhi(iOP)));
                    for iF = 1:3
                        % Plot arches
                        fill(...
                            OP_meanX(iOP) + x(iF,:),...
                            OP_meanY(iOP) + y(iF,:),EnKF.Colors(1,:),...
                            'FaceAlpha',0.1,'EdgeColor','none')
                        
                    end
                    % plot a line from OP to center of Arch
                    plot([OP_meanX(iOP),...
                        OP_meanX(iOP) + cos(phiW)*dataU(iOP)],...
                        [OP_meanY(iOP),...
                        OP_meanY(iOP) + sin(phiW)*dataU(iOP)],...
                        'Color',EnKF.Colors(5,:),'LineWidth',0.5)
                    
                end
                hold off
            elseif varU_Phi(1)
                % Only Variance in U
                
            elseif varU_Phi(2)
                % Only Variance in phi
                hold on
                for iOP = 1:length(dataU)
                    phiW = angSOWFA2world(dataPhi(iOP));
                    [x,y] = getPizzaSlice(dataU(iOP),...
                        phiW,deg2rad(stdPhi(iOP)));
                    for iF = 1:3
                        % Plot arches
                        fill(...
                            OP_meanX(iOP) + x(iF,:),...
                            OP_meanY(iOP) + y(iF,:),EnKF.Colors(1,:),...
                            'FaceAlpha',0.1,'EdgeColor','none')
                        
                    end
                    % plot a line from OP to center of Arch
                    plot([OP_meanX(iOP),...
                        OP_meanX(iOP) + cos(phiW)*dataU(iOP)],...
                        [OP_meanY(iOP),...
                        OP_meanY(iOP) + sin(phiW)*dataU(iOP)],...
                        'Color',EnKF.Colors(5,:),'LineWidth',0.5)
                    
                end
                hold off
            else
                return
            end
            grid on
            axis equal
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
        end
        
        function obj = plotK_WeightedFlowAndSOWFA(obj, EnKF, T, Vis, Sim, paramFLORIS, Wind,windDir)
            %% Get figure
            if isempty(obj.Figs{32});obj.Figs{32} = figure; end
            try
                figure(obj.Figs{32})
            catch
                obj.Figs{32} = figure;
            end
            clf;
            gap_HW    = [.06 .08];
            marg_h_LU = [.06 .06];
            marg_w_LR = [.06 .04];
            ha = tight_subplot(2,3,gap_HW,marg_h_LU,marg_w_LR);
            f = gcf;
            f.Position(3:4) = [757 592];

            %% Create 2x2 layout
            % 11 -> FLORIDyn estimate based on get tmpWF at SOWFA loactions
            % 12 -> SOWFA
            % 21 -> FLORIDyn variance
            % 22 -> Weights



            %% Get SOWFA data
            path2VTK = [Vis.FlowField.Error.ValidationPath ...
                num2str(Sim.EndTime) filesep Vis.FlowField.Error.SnapName];
            filename = [num2str(Sim.EndTime) 'FlowError'];
            
            if strcmp(path2VTK(end-4:end),'.vtk')
                [~,cellCenters,cellData] = importVTK(path2VTK);
            else
                [~,cellCenters,cellData] = importVTK([path2VTK '.vtk']);
            end
            
            Validation_UmeanAbs = sqrt(sum(cellData.^2,2));
            H = mean(cellCenters(:,3));
            F = scatteredInterpolant(cellCenters(:,1),cellCenters(:,2),...
                Validation_UmeanAbs);


            %% Combine Ensemble states and create estimate
            T.States_OP = mean(...
                reshape(EnKF.States_OP,...
                [T.nOP*T.nT,length(T.Names_OP),EnKF.nE])...
                ,3);
            
            T.States_WF(:,1) = mean(EnKF.States.Vel,2);
            T.States_WF(:,2) = mean(EnKF.States.Dir,2);

            % Set up wind field
            nM = 3;
            fieldLims = Vis.FlowField.Lims;
            fieldRes  = Vis.FlowField.Res;
            xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
            yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);

            [X,Y] = meshgrid(xAx,yAx);
            clear fieldRes xAx yAx
            
            Z = getMeasurementsPar2(X, Y, T.posNac(1,3), T, paramFLORIS, Wind);
            
            W_mat = getWeights(X(:), Y(:), T, Sim.Dyn.Vel, Sim.TimeStep);
            V = reshape(W_mat./sum(W_mat,2) * T.States_WF(:,1),size(X));

            %% Weights
            

            %% Variance
            % run estimate + 1 sigma and calc difference
%             T.States_WF(:,1) = mean(EnKF.States.Vel,2) + ...
%                 std(EnKF.States.Vel,[],2);
%             Z_V = getMeasurementsPar2(X, Y, T.posNac(1,3), T, paramFLORIS, Wind);

            %% Plot wind speed
            axes(ha(1))
            contourf(X,Y,V .* Z(:,:,1),21,'LineColor','none')
            hold on
            % Plot Rotors
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                %yaw = deg2rad(T.States_T(T.StartI(i_T),2));
                yaw = angSOWFA2world(windDir- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'k','LineWidth',3);
            end

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[1 1 1])

            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Wind speed (ms^{-1})';
            colormap(ha(1),viridis(16))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            caxis([2,10])
            title('Velocity based on EnKF')

            %% Plot variance
            axes(ha(4))
            contourf(X,Y,reshape(W_mat./sum(W_mat,2)*...
                std(EnKF.States.Vel,[],2),size(X)),21,'LineColor','none')
            hold on
            % Plot Rotors
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                %yaw = deg2rad(T.States_T(T.StartI(i_T),2));
                yaw = angSOWFA2world(windDir- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'k','LineWidth',3);
            end

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[1 1 1])
            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Wind speed std (ms^{-1})';
            colormap(ha(4),inferno(16))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            caxis([0.1,.3])
            title('Velocity std based on EnKF')
            
            %%
            axes(ha(2))
            dataDir = reshape(W_mat./sum(W_mat,2)*...
                T.States_WF(:,2),size(X));
            contourf(X,Y,dataDir,21,'LineColor','none')
            hold on
            % Plot Rotors
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                %yaw = deg2rad(T.States_T(T.StartI(i_T),2));
                yaw = angSOWFA2world(windDir- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'k','LineWidth',3);
            end

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[1 1 1])

            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Wind direction (deg)';
            colormap(ha(2),RdBu(17))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')

            mData = mean(dataDir,"all");
            difMM = max(abs(dataDir),[],'all') - mData;
            caxis(windDir + [-0.5,.5])
            title('Direction based on EnKF')

            %% Wind dir variance
            axes(ha(5))
            contourf(X,Y,reshape(W_mat./sum(W_mat,2)*...
                std(EnKF.States.Dir,[],2),size(X)),21,'LineColor','none')
            hold on
            % Plot Rotors
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                %yaw = deg2rad(T.States_T(T.StartI(i_T),2));
                yaw = angSOWFA2world(windDir- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'k','LineWidth',3);
            end

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[1 1 1])
            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Wind direction std (deg)';
            colormap(ha(5),inferno(16))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            caxis([0.3,0.9])
            title('Direction std based on EnKF')
            %% Plot SOWFA
            axes(ha(3))
            contourf(X,Y,F(X,Y),21,'LineColor','none')
            hold on
            % Plot Rotors
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                %yaw = deg2rad(T.States_T(T.StartI(i_T),2));
                yaw = angSOWFA2world(windDir- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'k','LineWidth',3);
            end

            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Wind speed (ms^{-1})';
            colormap(ha(3),viridis(16))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            caxis([2,10])
            title('SOWFA flow field')

            %% Plot SOWFA w FLORIDyn
            axes(ha(6))
            contourf(X,Y,F(X,Y),21,'LineColor','none')
            hold on
            % Plot Rotors
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                %yaw = deg2rad(T.States_T(T.StartI(i_T),2));
                yaw = angSOWFA2world(windDir- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'k','LineWidth',3);
            end

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[1 1 1])
            contour(X,Y,V .* Z(:,:,1),10,'LineColor',[1,1,1])

            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Wind speed (ms^{-1})';
            colormap(ha(6),viridis(16))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            caxis([2,10])
            title('SOWFA with FLORIDyn contour')

            %% Store
            folderName = [Sim.PathToSim filesep 'Results' ...
                        filesep 'EnKF_weighted_FlowField' filesep 'Comparison'];
            if not(isfolder(folderName))
                mkdir(folderName)
            end
            saveas(f,[folderName filesep filename '.png'])
        end

        function obj = plotK_WeightedFlowField(obj, EnKF, T, Vis, Sim, iState, iS, varargin)
            if isempty(obj.Figs{12}); obj.Figs{12} = figure; end
            try
                figure(obj.Figs{12})
            catch
                obj.Figs{12} = figure;
            end

            clf;
            gap_HW    = [.08 .06];
            marg_h_LU = [.06 .06];
            marg_w_LR = [.06 .04];
            ha = tight_subplot(1,4,gap_HW,marg_h_LU,marg_w_LR);
            f = gcf;
            f.Position = [150 190 1155*4/3 390];
            
            x = reshape(EnKF.States_OP(:,1:EnKF.nStatesOP:end),[],1);
            y = reshape(EnKF.States_OP(:,2:EnKF.nStatesOP:end),[],1);
            
            fieldLims = Vis.FlowField.Lims;
            
            inBounds = x<fieldLims(2,1)+100;
            inBounds = and(y<fieldLims(2,2)+100,inBounds);
            inBounds = and(x>fieldLims(1,1)-100,inBounds);
            inBounds = and(y>fieldLims(1,2)-100,inBounds);
            x = x(inBounds);
            y = y(inBounds);
            
            fieldRes  = Vis.FlowField.Res;
            xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
            yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);
            
            [X,Y] = meshgrid(xAx,yAx);
            if nargin > 7
                normVals = true;
            else
                normVals = false;
            end
            
            switch iState
                case 1
                    val = mean(EnKF.States.Vel,2);
                    if normVals
                        val = varargin{1} - val;
                        borderCmap = [-2,2];
                        cLabel = ['Wind speed difference to '...
                            num2str(varargin{1}) ' ms^{-1}'];

                        valVar = var(EnKF.States.Vel,[],2)./varargin{1};
                        borderCmapVar = [0,0.1];
                        cLabelVar = 'Norm. Variance (-)';
                    else
                        borderCmap = [7,10];
                        cLabel = 'Wind speed (ms^{-1})';

                        valVar = var(EnKF.States.Vel,[],2);
                        borderCmapVar = [0,3];
                        cLabelVar = 'Variance (ms^{-1})';
                    end
                    
                    folderName = [Sim.PathToSim filesep 'Results' ...
                        filesep 'EnKF_weighted' filesep 'WindSpeed'];

                    Dyn = Sim.Dyn.Vel;
                case 2
                    val = mean(EnKF.States.Dir,2);
                    valVar = var(EnKF.States.Dir,[],2);
                    if normVals
                        val = varargin{1} - val;
                        valVar = var(EnKF.States.Dir,[],2)./360;
                        borderCmap = [-5,5];
                        cLabel = ['Wind direction difference to '...
                            num2str(varargin{1}) ' deg'];
                        borderCmapVar = [0,0.1];
                        cLabelVar = 'Variance normed by 360 deg';
                    else
                        borderCmap = [180,270];
                        cLabel = 'Wind direction (deg)';
                        borderCmapVar = [0,50];
                        cLabelVar = 'Variance (deg)';
                    end

                    folderName = [Sim.PathToSim filesep 'Results' ...
                        filesep 'EnKF_weighted' filesep 'WindDirection'];

                    Dyn = Sim.Dyn.Dir;
                case 3
                    val = mean(EnKF.States.TI,2);
                    val = val(inBounds);
                    borderCmap = [0,0.3];
                    cLabel = 'Turbulence intensity (%)';

                    folderName = [Sim.PathToSim filesep 'Results' ...
                        filesep 'EnKF_weighted' filesep 'TI'];
            end

            %% get weights
            W_mat = getWeights(X(:), Y(:), T, Dyn, Sim.TimeStep);
            
            %% Plot
            %W = reshape(sum(W_mat,2),size(X));
            W = reshape(max(W_mat,[],2),size(X));
            V = reshape(W_mat./sum(W_mat,2) * val,size(X));
            Vvar = reshape(W_mat./sum(W_mat,2) * valVar,size(X));

            axes(ha(1))
            contourf(X,Y,V,21,'EdgeColor','none');
            if normVals
                colormap(ha(1),flipud(RdBu(21)))
            else
                colormap(ha(1),viridis(21))
            end
            c = colorbar('southoutside');
            c.Label.String = cLabel;
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[361,361],'k','LineWidth',2);
            end
            scatter3(x,y,ones(size(x))*361,1,zeros(1,3),'filled')
            hold off
            title('Weighted mean data')
            caxis(borderCmap)
            
            axes(ha(2))
            contourf(X,Y,Vvar,21,'EdgeColor','none');
            colormap(ha(2),inferno(21))
            c = colorbar('southoutside');
            c.Label.String = cLabelVar;
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[361,361],'w','LineWidth',2);
            end
            scatter3(x,y,ones(size(x))*361,1,ones(1,3),'filled')
            hold off
            title('Weighted variance data')
            caxis(borderCmapVar)

            % Plot weighting function
            axes(ha(3))
            contourf(X,Y,W,10,'EdgeColor','none');
            colormap(ha(3),inferno(10))
            c = colorbar('southoutside');
            c.Label.String = 'Weight (-)';
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'w','LineWidth',2);
            end
            scatter3(x,y,ones(size(x))*361,1,ones(1,3),'filled')
            hold off
            caxis([0,1])
            title('Maximum weight')
            
            % Plot values with transparency 
            axes(ha(4))
            s = surf(X,Y,V,'EdgeColor','none');
            hold on
            alpha(s,tanh((rescale(W)-.5)*5)*.5+.5)
            contour3(X,Y,W+360,5,'k')
            hold off
            if normVals
                colormap(ha(4),flipud(RdBu(21)))
            else
                colormap(ha(4),viridis(21))
            end
            c = colorbar('southoutside');
            c.Label.String = cLabel;
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[361,361],'k','LineWidth',2);
            end
            hold off
            title('Trusted data')
            caxis(borderCmap)
            
            if not(isfolder(folderName))
                mkdir(folderName)
            end
            saveas(f,[folderName filesep num2str(iS) '.png'])
        end

        function obj = plotK_GaussianFlowField(obj, EnKF, T, Vis, Sim, iState, iS)
            %% Plot the entire flow field weighted by a Gauss function, 
            % the function provides weights to all OPs, calculated by the
            % distance to the Grid point of interest.
            % First figure shows all Grid points and their raw data
            % Second figure shows the sum of weights that went into the
            % grid point - the higher, the more OPs contributed, the more
            % certainty we have in the value
            % Third figure shows the raw data but with transpareny based on
            % the weights from the second figure
            if isempty(obj.Figs{12}); obj.Figs{12} = figure; end
            try
                figure(obj.Figs{12})
            catch
                obj.Figs{12} = figure;
            end
            
            if iS == 1
                setLims = true;
            else
                setLims = false;
            end
            
            clf;
            gap_HW    = [.08 .06];
            marg_h_LU = [.06 .06];
            marg_w_LR = [.06 .04];
            ha = tight_subplot(1,3,gap_HW,marg_h_LU,marg_w_LR);
            f = gcf;
            f.Position = [520 190 1155 390];
            
            x = reshape(EnKF.States_OP(:,1:EnKF.nStatesOP:end),[],1);
            y = reshape(EnKF.States_OP(:,2:EnKF.nStatesOP:end),[],1);
            
            fieldLims = Vis.FlowField.Lims;
            
            inBounds = x<fieldLims(2,1)+100;
            inBounds = and(y<fieldLims(2,2)+100,inBounds);
            inBounds = and(x>fieldLims(1,1)-100,inBounds);
            inBounds = and(y>fieldLims(1,2)-100,inBounds);
            x = x(inBounds);
            y = y(inBounds);
            
            fieldRes  = Vis.FlowField.Res;
            xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
            yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);
            
            [X,Y] = meshgrid(xAx,yAx);
            V = zeros(size(X));
            W = zeros(size(X));
            
            switch iState
                case 1
                    val = EnKF.States.Vel(:);
                    val = val(inBounds);
                    borderCmap = 10;
                    cLabel = 'Wind speed (m/s)';
                case 2
                    val = EnKF.States.Dir(:);
                    val = val(inBounds);
                    borderCmap = 20;
                    cLabel = 'Wind direction (deg)';
                case 3
                    val = EnKF.States.TI(:);
                    val = val(inBounds);
                    borderCmap = 0.4;
                    cLabel = 'Turbulence intensity (%)';
            end
            
            
            parfor iGP = 1:length(X(:))
                [V(iGP), W(iGP)] = scatteredGaussian(...
                    x,y,val,X(iGP),Y(iGP),T.D(1)/2);
                
            end
            % TMP FIX
            threshW = 1.5;
            W(W>threshW) = threshW;
            
            % Set adaptive limits
            maxV = max(V,[],'all');
            minV = min(V,[],'all');
            
            if setLims
                obj.tmpLims = [minV maxV] + [-.1 .1]*max((maxV-minV),borderCmap);
            else
                obj.tmpLims = obj.tmpLims + ...
                    0.5*(([minV maxV] + [-.1 .1]*max((maxV-minV),borderCmap)) - obj.tmpLims);
            end
            
            % Plot Values
            axes(ha(1))
            surf(X,Y,V,'EdgeColor','none');
            colormap(ha(1),viridis(1000))
            c = colorbar('southoutside');
            c.Label.String = cLabel;
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[361,361],'w','LineWidth',2);
            end
            scatter3(x,y,ones(size(x))*361,1,ones(1,3),'filled')
            hold off
            title('Gaussian weighted data')
            caxis(obj.tmpLims)
            
            % Plot weighting function
            axes(ha(2))
            surf(X,Y,W,'EdgeColor','none');
            colormap(ha(2),inferno(1000))
            c = colorbar('southoutside');
            c.Label.String = 'Weight (-)';
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'w','LineWidth',2);
            end
            scatter3(x,y,ones(size(x))*361,1,ones(1,3),'filled')
            hold off
            %caxis([0,1])
            title('Collective weights')
            
            % Plot values with transparency 
            axes(ha(3))
            s = surf(X,Y,V,'EdgeColor','none');
            alpha(s,tanh((rescale(W)-.5)*5)*.5+.5)
            colormap(ha(3),viridis(1000))
            c = colorbar('southoutside');
            c.Label.String = cLabel;
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[361,361],'k','LineWidth',2);
            end
            hold off
            title('Trusted data')
            caxis(obj.tmpLims)
            
            if not(isfolder([Sim.PathToSim filesep 'Results' filesep 'EnKF']))
                mkdir([Sim.PathToSim filesep 'Results' filesep 'EnKF'])
            end
            saveas(f,[Sim.PathToSim 'Results' filesep 'EnKF' filesep num2str(iS) '.png'])
        end
        
        function [obj, err] = plotK_FlowField(obj, EnKF, T, Wind,Sim,Vis,paramFLORIDyn,paramFLORIS,iState,trueDir)
            if isempty(obj.Figs{25+iState});obj.Figs{25+iState} = figure; end
            try
                figure(obj.Figs{25+iState})
            catch
                obj.Figs{25+iState} = figure;
            end
            clf;
            gap_HW    = [.08 .06];
            marg_h_LU = [.06 .06];
            marg_w_LR = [.06 .04];
            ha = tight_subplot(1,3,gap_HW,marg_h_LU,marg_w_LR);
            f = gcf;
            f.Position(3:4) = [1155 390];

            %% Combine EnKF states to mean model
            T.States_OP = mean(...
                reshape(EnKF.States_OP,...
                [T.nOP*T.nT,length(T.Names_OP),EnKF.nE])...
                ,3);
            
            if EnKF.Vel.Correct
                T.States_WF(:,1) = mean(EnKF.States.Vel,2);
                u_std = std(EnKF.States.Vel,[],2);
            end
            if EnKF.Dir.Correct
                T.States_WF(:,2) = mean(EnKF.States.Dir,2);
                dir_std = std(EnKF.States.Dir,[],2);
            end
            if EnKF.TI.Correct
                T.States_WF(:,3) = mean(EnKF.States.TI,2);
            end
            %% =========== Get reduction
            nM = 3;
            fieldLims = Vis.FlowField.Lims;
            fieldRes  = Vis.FlowField.Res;
            xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
            yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);

            [X,Y] = meshgrid(xAx,yAx);
            clear fieldRes xAx yAx

            %% =========== Get wind field data

%             if Vis.FlowField.Plot.parallel
                Z = getMeasurementsPar2(X, Y, T.posNac(1,3), T, paramFLORIS, Wind);
%             else
%                 Z = getMeasurements(X,Y,nM,T.posNac(1,3),T,paramFLORIS,Wind,Vis);
%             end
            
            %% =========== get weights
            W_mat = getWeights(X(:), Y(:), T, Sim.Dyn.Vel, Sim.TimeStep);
            V = reshape(W_mat./sum(W_mat,2) * T.States_WF(:,1),size(X));

            %% =========== Load vtk data ================================ %
            path2VTK = [Vis.FlowField.Error.ValidationPath ...
                num2str(Sim.EndTime) filesep Vis.FlowField.Error.SnapName];
            filename = [num2str(Sim.EndTime) 'FlowError'];
            
            if strcmp(path2VTK(end-4:end),'.vtk')
                [~,cellCenters,cellData] = importVTK(path2VTK);
            else
                [~,cellCenters,cellData] = importVTK([path2VTK '.vtk']);
            end
            
            Validation_UmeanAbs = sqrt(sum(cellData.^2,2));
            H = mean(cellCenters(:,3));
            F = scatteredInterpolant(cellCenters(:,1),cellCenters(:,2),...
                Validation_UmeanAbs);

            %% Plot combined
            axes(ha(1))
            contourf(X,Y,V .* Z(:,:,1),[-20,(2:.5:10)],'LineColor','none')
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

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[1 1 1])

            scatter(T.States_OP(1:10:end,1),T.States_OP(1:10:end,2),15,...
                'filled','MarkerFaceColor',[1 1 1])

            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Wind speed (ms^{-1})';
            colormap(ha(1),viridis(16))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            caxis([2,10])
            title('EnKF based FLORIDyn state')
            
            %% Plot error
            axes(ha(2))
            contourf(X,Y,F(X,Y) - V .* Z(:,:,1),...
                [-20,linspace(-3.1,3.1,31)], 'LineColor','none')
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

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[0 0 0])

            scatter(T.States_OP(1:10:end,1),T.States_OP(1:10:end,2),15,...
                'filled','MarkerFaceColor',[0 0 0])

            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Diff. wind speed (ms^{-1})';
            colormap(ha(2),RdBu(31))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            caxis([-3.1,3.1])
            title('Error between FLORIDyn and SOWFA')
            %% Plot original
            axes(ha(3))
            contourf(X,Y,F(X,Y),[-20,(2:.5:10)],'LineColor','none')
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

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[1 1 1])

            scatter(T.States_OP(1:10:end,1),T.States_OP(1:10:end,2),15,...
                'filled','MarkerFaceColor',[1 1 1])

            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Wind speed (ms^{-1})';
            colormap(ha(3),viridis(16))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            caxis([2,10])
            title('SOWFA flow field')
            %% Store
%             folderName = [Sim.PathToSim filesep 'Results' ...
%                         filesep 'EnKF_weighted' filesep 'Validation'];
            folderName = [Sim.PathToResults 'Validation'];
            if not(isfolder(folderName))
                mkdir(folderName)
            end
            saveas(f,[folderName filesep filename '.png'])
            
            %% Plot FLORIDyn standalone
            f = figure(2929);
            f.Position(3:4) = [500 500];
            contourf(X,Y,V .* Z(:,:,1),[-20,(2:.5:10)],'LineColor','none')
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

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[1 1 1])

            scatter(T.States_OP(1:10:end,1),T.States_OP(1:10:end,2),15,...
                'filled','MarkerFaceColor',[1 1 1])

            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Wind speed (ms^{-1})';
            colormap(viridis(16))
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            caxis([2,10])
            title(['EnKF based FLORIDyn state, t=' num2str(Sim.EndTime) ' s'])
            %% Store
%             folderName = [Sim.PathToSim filesep 'Results' ...
%                         filesep 'Movie' filesep 'EnKF'];
            folderName = [Sim.PathToResults 'Snapshot'];
            if not(isfolder(folderName))
                mkdir(folderName)
            end
            saveas(f,[folderName filesep  'Mov_t' num2str(Sim.EndTime) '.png'])

            %% Error calculation
            f = figure(292929);
            f.Position(3:4) = [500 500];
            Std = reshape(W_mat./sum(W_mat,2) * u_std,size(X));
            Err = F(X,Y) - V .* Z(:,:,1);
            ErrN = Err./Std;
            contourf(X,Y,ErrN,[-1000,(-4:4)],'LineColor','none')
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

            scatter(T.States_OP(:,1),T.States_OP(:,2),5,...
                'filled','MarkerFaceColor',[1 1 1]*0)

            scatter(T.States_OP(1:10:end,1),T.States_OP(1:10:end,2),15,...
                'filled','MarkerFaceColor',[1 1 1]*0)

            axis equal
            xlim([fieldLims(1,1),fieldLims(2,1)])
            ylim([fieldLims(1,2),fieldLims(2,2)])

            hold off

            c = colorbar('southoutside');
            c.Label.String ='Error norm. by std.(-)';
            colormap(RdBu(9))
            caxis([-4.5,4.5])
%             xlabel('West-East (m)')
%             ylabel('South-North (m)')
            title(['Norma. error SOWFA to EnKF state, t=' num2str(Sim.EndTime) ' s'])

            %% Store
%             folderName = [Sim.PathToSim filesep 'Results' ...
%                         filesep 'NormErr' filesep 'EnKF'];
            folderName = [Sim.PathToResults 'NormErr'];
            if not(isfolder(folderName))
                mkdir(folderName)
            end
            saveas(f,[folderName filesep  'NormErr_t' num2str(Sim.EndTime) '.png'])
            
            %% Calculate OP positional variance
            stdOPs_3d = std(...
                reshape(EnKF.States_OP,...
                [T.nOP*T.nT,length(T.Names_OP),EnKF.nE]),[],3);

            stdOPs = sqrt(sum(stdOPs_3d(:,1:2).^2,2));
            inBounds = and(T.States_OP(:,1)<Vis.FlowField.Lims(2,1),...
                T.States_OP(:,2)<Vis.FlowField.Lims(2,2));
            
            %% Caluclate wind speed error
            V = reshape(W_mat./sum(W_mat,2) * T.States_WF(:,2),size(X));
            ErrD = trueDir - V;
            ErrDN = ErrD./reshape(W_mat./sum(W_mat,2) * dir_std,size(X));
            
            %% Error output
            % Error stats
            % 1:  Average RMSE of the wind speed across the field
            % 2:  Average Error of the wind speed across the field
            % 3:  Std of the Error of the wind speed across the field
            % 4:  Percentage of wind speed errors within 1 std
            % 5:  Percentage of wind speed errors within 2 std
            % 6:  Percentage of wind speed errors within 3 std
            % ==== OP Pos. ====
            % 7:  Average OP position std
            % 8:  Average OP position std within wind farm bounds
            % 9:  Average OP position std outside wind farm bounds
            % ==== Wind Dir ====
            % 10: Average RMSE of the wind direction across the field
            % 11: Average Error of the wind direction across the field
            % 12: Std of the Error of the wind direction across the field
            % 13: Percentage of wind direction errors within 1 std
            % 14: Percentage of wind direction errors within 2 std
            % 15: Percentage of wind direction errors within 3 std
            err = [...
                sum(sqrt(Err(:).^2))/length(Err(:)),...
                mean(Err(:)),...
                std(Err(:)), ...
                sum(abs(ErrN(:))<1)/length(Err(:)), ...
                sum(abs(ErrN(:))<2)/length(Err(:)), ...
                sum(abs(ErrN(:))<3)/length(Err(:)), ...
                mean(stdOPs),...
                mean(stdOPs(inBounds)),...
                mean(stdOPs(~inBounds)),...
                sum(sqrt(ErrD(:).^2))/length(ErrD(:)),...
                mean(ErrD(:)),...
                std(ErrD(:)), ...
                sum(abs(ErrDN(:))<1)/length(ErrD(:)), ...
                sum(abs(ErrDN(:))<2)/length(ErrD(:)), ...
                sum(abs(ErrDN(:))<3)/length(ErrD(:))];
        end
        
        function obj = plotK_Correction(obj, EnKF, T, Vis, K, diff, iE, iS, iState)
            %% Plot how much which turbine corrects which state
            if isempty(obj.Figs{13+iState});obj.Figs{13+iState} = figure; end
            try
                figure(obj.Figs{13+iState})
            catch
                obj.Figs{13+iState} = figure;
            end
            
            if and(iS == 1, iE == 1)
                clf;
                l = getLayout(T.nT);
                gap_HW    = [.01 .01];
                marg_h_LU = [.06 .06];
                marg_w_LR = [.06 .06];
                obj.AxesHandles{13+iState} = ...
                    tight_subplot(l(1),l(2),gap_HW,marg_h_LU,marg_w_LR);
            end
            if iE>1; return; end % <= for now
            
            ha = obj.AxesHandles{13+iState};
            f.Position = [520 190 830 650];
            
            for iT = 1:T.nT
                axes(ha(iT))
                hold on
                imagesc([1 T.nOP*T.nT],[iS iS],K(:,iT)')%(K(:,iT)*diff(iT))')
                hold off
                xlim([1 T.nOP*T.nT])
                ylim([1 EnKF.Sim.Sections])
                colormap(flipud(RdBu))
                caxis([-1,1])
            end
            
            
        end
        
        function obj = plotK_Correction_Quiver(obj, EnKF, T, Vis, K, diff, iE, iS, iState, Sim)
            %% Plot how much which turbine corrects which state
            if isempty(obj.Figs{22+iState});obj.Figs{22+iState} = figure; end
            try
                figure(obj.Figs{22+iState})
            catch
                obj.Figs{22+iState} = figure;
            end
            
            if iE>1; return; end % <= for now
            f = gcf;
            if iS == 1
                clf;
            end
            
            f.Position = [520 190 1.3552e+03 650];
            cMax = max(abs(K.*repmat(diff',T.nT*T.nOP,1)),[],'all');
            cmap = flipud(RdBu(101));
            subplot(1,2,2)
            meanCorr = mean(K.*repmat(diff',T.nT*T.nOP,1),'all');
            stdCorr = std(K.*repmat(diff',T.nT*T.nOP,1),[],'all');
            
            hold on
            scatter(iS,meanCorr,20,EnKF.Colors(5,:),'filled')
            plot([iS,iS], [-1,1]*stdCorr + meanCorr,'_-','Color', EnKF.Colors(5,:))
            grid on
            xlim([0,EnKF.Sim.Sections+1])
            
            ha = subplot(1,2,1);
            cla(ha)
            for iT = 1:T.nT
                xT = repmat(EnKF.States_OP(T.StartI(iT),(iE-1)*EnKF.nStatesOP+1),T.nOP*T.nT,1);
                yT = repmat(EnKF.States_OP(T.StartI(iT),(iE-1)*EnKF.nStatesOP+2),T.nOP*T.nT,1);
                xD = EnKF.States_OP(:,(iE-1)*EnKF.nStatesOP+1) - xT;
                yD = EnKF.States_OP(:,(iE-1)*EnKF.nStatesOP+2) - yT;
                abThresh = abs(K(:,iT)*diff(iT)) > 0.1;
                xT = xT(abThresh);
                yT = yT(abThresh);
                xD = xD(abThresh);
                yD = yD(abThresh);
                idx = round(...
                    max(min(K(abThresh,iT)*diff(iT)/cMax,1),-1) * 50) + 51;
                if iT == 1; hold on; end
                for iQ = 1:length(xT)
                    quiver(xT(iQ),yT(iQ),xD(iQ),yD(iQ),...
                        'AutoScale','off','Color',cmap(idx(iQ),:));
                end
            end
            hold off
            grid on
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            title('Correction')
            
            if not(isfolder([Sim.PathToSim filesep 'Results' filesep 'EnKF_QCorrection']))
                mkdir([Sim.PathToSim '\Results' filesep 'EnKF_QCorrection'])
            end
            saveas(f,[Sim.PathToSim 'Results' filesep 'EnKF_QCorrection' filesep num2str(iS) '.png'])
        end
        
        function obj = plotK_GaussianFlowFieldError(obj, EnKF, T, Vis, Sim, iS, paramFLORIS, Wind)
            %% Plot the entire flow field error weighted by a Gauss 
            % function, 
            if ~obj.Settings.FlowFieldError.Online
                dirSnap = dir(Vis.FlowField.Error.ValidationPath);
                Vis.FlowField.Error.Steps = zeros(length(dirSnap)-2,1);
                for i = 3:length(dirSnap)
                    try
                        Vis.FlowField.Error.Steps(i-2) = ...
                            str2num(dirSnap(i).name);
                    catch
                        Vis.FlowField.Error.Steps(i-2) = -1;
                    end
                end
                if sum(Vis.FlowField.Error.Steps==Sim.EndTime)==0;return;end
            end
            
            
            % Can only compare velocity
            if ~EnKF.Vel.Correct; return; end
            
            % Init figure
            if isempty(obj.Figs{22}); obj.Figs{22} = figure; end
            try
                figure(obj.Figs{22})
            catch
                obj.Figs{22} = figure;
            end
            
            if iS == 1
                setLims = true;
            else
                setLims = false;
            end
            % ============ Load vtk data ================================ %
            path2VTK = [Vis.FlowField.Error.ValidationPath ...
                num2str(Sim.EndTime) '\' Vis.FlowField.Error.SnapName];
            filename = [num2str(Sim.EndTime) 'FlowError'];
            
            if strcmp(path2VTK(end-4:end),'.vtk')
                [~,cellCenters,cellData] = importVTK(path2VTK);
            else
                [~,cellCenters,cellData] = importVTK([path2VTK '.vtk']);
            end
            
            Validation_UmeanAbs = sqrt(sum(cellData.^2,2));
            H = mean(cellCenters(:,3));
            F = scatteredInterpolant(cellCenters(:,1),cellCenters(:,2),...
                Validation_UmeanAbs);
            
            % ============ Prep figure ================================== %
            clf;
            gap_HW    = [.08 .06];
            marg_h_LU = [.06 .06];
            marg_w_LR = [.06 .04];
            ha = tight_subplot(1,3,gap_HW,marg_h_LU,marg_w_LR);
            f = gcf;
            f.Position = [520 190 1155 390];
            
            % ============ Make Mesh ==================================== %
%             x = reshape(EnKF.States_OP(:,1:EnKF.nStatesOP:end),[],1);
%             y = reshape(EnKF.States_OP(:,2:EnKF.nStatesOP:end),[],1);
%             
            fieldLims = Vis.FlowField.Lims;
%             
%             inBounds = x<fieldLims(2,1)+100;
%             inBounds = and(y<fieldLims(2,2)+100,inBounds);
%             inBounds = and(x>fieldLims(1,1)-100,inBounds);
%             inBounds = and(y>fieldLims(1,2)-100,inBounds);
%             x = x(inBounds);
%             y = y(inBounds);
            
            fieldRes  = Vis.FlowField.Res;
            xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
            yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);
            
            [X,Y] = meshgrid(xAx,yAx);
            V = zeros(size(X));
            W = zeros(size(X));
            
            % ============ Calc OP data ================================= %
%             val = EnKF.States.Vel(:);
%             val = val(inBounds);
            cLabel = 'Wind speed error (m/s)';
            
            nS = length(T.Names_OP);
            for iStates = 1:nS
                T.States_OP(:,iStates) = ...
                    mean(EnKF.States_OP(:,iStates:nS:end),2);
            end
            T.States_WF(:,1) = EnKF.States.Vel_mean;
            if EnKF.Dir.Correct
                T.States_WF(:,2) = EnKF.States.Dir_mean;
            end
            if EnKF.TI.Correct
                T.States_WF(:,3) = EnKF.States.TI_mean;
            end
            
            V = getMeasurementsPar2(X, Y, H, T, paramFLORIS, Wind);
            
            val = std(EnKF.States.Vel,[],2);
            val = 1-rescale(val);
            
            TD = T.D(1); 
            x = T.States_OP(:,1);
            y = T.States_OP(:,2);
            GaussSOWFA = zeros(size(X));
            parfor iGP = 1:length(X(:))
                [W(iGP), ~] = scatteredGaussian(...
                    x,y,val,...
                    X(iGP),Y(iGP),TD/2);
                [GaussSOWFA(iGP), ~] = scatteredGaussian(...
                    cellCenters(:,1),cellCenters(:,2),...
                    Validation_UmeanAbs,...
                    X(iGP),Y(iGP),TD/4);
                
            end
            % ============ Calc difference ============================== %
            %V = F(X,Y) - V;
            V = GaussSOWFA - V;
            % Set adaptive limits
            maxV = max(abs(V),[],'all');
            minV = -maxV;
            maxV = 4; minV = -4; % <=================== Constant overwrite
            % Plot Values
            axes(ha(1))
            contourf(X,Y,V,linspace(minV,maxV,11),'EdgeColor','none')
            %surf(X,Y,V,'EdgeColor','none');
            colormap(ha(1),flipud(RdBu(11)))
            c = colorbar('southoutside');
            c.Label.String = cLabel;
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[361,361],'w','LineWidth',2);
            end
            scatter3(x,y,ones(size(x))*361,1,zeros(1,3),'filled')
            hold off
            title('Flow field based on mean states')
            caxis([minV maxV])
            
            % Plot weighting function
            axes(ha(2))
            surf(X,Y,rescale(W),'EdgeColor','none');
            colormap(ha(2),inferno(1000))
            c = colorbar('southoutside');
            c.Label.String = 'Weight (-)';
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),2));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[20,20],'w','LineWidth',2);
            end
            scatter3(x,y,ones(size(x))*361,1,ones(1,3),'filled')
            hold off
            caxis([0,1])
            title('Certainty based on std')
            
            % Plot values with transparency 
            axes(ha(3))
            s = surf(X,Y,V,'EdgeColor','none');
            alpha(s,rescale(W))
            colormap(ha(3),flipud(RdBu(1000)))
            c = colorbar('southoutside');
            c.Label.String = cLabel;
            grid on
            axis equal
            view(2)
            xlim(Vis.FlowField.Lims(:,1)')
            ylim(Vis.FlowField.Lims(:,2)')
            hold on
            for i_T = 1:length(T.D)
                % Get start and end of the turbine rotor
                yaw = angSOWFA2world(T.States_WF(T.StartI(i_T),2)- ...
                    T.States_T(T.StartI(i_T),1));
                rot_pos = ...
                    [cos(yaw), -sin(yaw);...
                    sin(yaw), cos(yaw)] * ...
                    [0,0;T.D(i_T)/2,-T.D(i_T)/2];
                rot_pos = rot_pos + repmat(T.posBase(i_T,1:2)',1,size(rot_pos,2));
                plot3(rot_pos(1,:),rot_pos(2,:),[361,361],'k','LineWidth',2);
            end
            hold off
            title('Trusted data')
            caxis([minV maxV])
            
            if not(isfolder([Sim.PathToSim filesep 'Results' filesep 'EnKF_FlowFieldError']))
                mkdir([Sim.PathToSim filesep 'Results' filesep 'EnKF_FlowFieldError'])
            end
            saveas(f,[Sim.PathToSim 'Results' filesep 'EnKF_FlowFieldError' filesep num2str(iS) '.png'])
        end
    end
end

function [x,y] = getArch(mU,sU,mD,sD)
%% Converts a mean Wind speed and variance along with mean wind dir and
% variance into a set of 5 overlapping arch segments
% Seg 1: +-3 sig wind speed and +-3 sig wind dir
phi = linspace(mD-3*sD,mD+3*sD,50);
x1 = [(mU+3*sU)*cos(phi),(mU-3*sU)*cos(fliplr(phi))];
y1 = [(mU+3*sU)*sin(phi),(mU-3*sU)*sin(fliplr(phi))];
% Seg 2: +-2 sig wind speed and +-2 sig wind dir
phi = linspace(mD-2*sD,mD+2*sD,50);
x2 = [(mU+2*sU)*cos(phi),(mU-2*sU)*cos(fliplr(phi))];
y2 = [(mU+2*sU)*sin(phi),(mU-2*sU)*sin(fliplr(phi))];
% Seg 3: +-1 sig wind speed and +-1 sig wind dir
phi = linspace(mD-1*sD,mD+1*sD,50);
x3 = [(mU+1*sU)*cos(phi),(mU-1*sU)*cos(fliplr(phi))];
y3 = [(mU+1*sU)*sin(phi),(mU-1*sU)*sin(fliplr(phi))];

x = [x1;x2;x3];
y = [y1;y2;y3];
x = [x,x(:,1)];
y = [y,y(:,1)];

% % Test code
% figure; hold on
% for iF = 1:3
%     fill(...
%         x(iF,:),...
%         y(iF,:),EnKF.Colors(1,:),...
%         'FaceAlpha',0.1,'EdgeColor','none')
%     % plot a line from OP to center of Arch
% end
% hold off
end

function [x,y] = getPizzaSlice(U,mD,sD)
%% Converts a mean Wind speed and variance along with mean wind dir and
% variance into a set of 5 overlapping arch segments
% Seg 1: +-3 sig wind speed and +-3 sig wind dir
phi = linspace(mD-3*sD,mD+3*sD,50);
x1 = [U*cos(phi),0];
y1 = [U*sin(phi),0];
% Seg 2: +-2 sig wind speed and +-2 sig wind dir
phi = linspace(mD-2*sD,mD+2*sD,50);
x2 = [U*cos(phi),0];
y2 = [U*sin(phi),0];
% Seg 3: +-1 sig wind speed and +-1 sig wind dir
phi = linspace(mD-1*sD,mD+1*sD,50);
x3 = [U*cos(phi),0];
y3 = [U*sin(phi),0];

x = [x1;x2;x3];
y = [y1;y2;y3];
x = [x,x(:,1)];
y = [y,y(:,1)];

% Test code
% figure; hold on
% for iF = 1:3
%     fill(...
%         x(iF,:),...
%         y(iF,:),EnKF.Colors(1,:),...
%         'FaceAlpha',0.1,'EdgeColor','none')
%     plot a line from OP to center of Arch
% end
% hold off
end

function [value, weight] = scatteredGaussian(x,y,v,xG,yG,sig)
%% Weighting function for the gaussian weighted plot
d = sqrt(sum((xG-x).^2+(yG-y).^2,2));
w = 1/(sig*sqrt(2*pi))*exp(-.5*(d/sig).^2);

weight = sum(w);
value = sum(w.*v)/weight;
end

function Z = getMeasurementsPar2(X,Y,zh,T,paramFLORIS,Wind)
% Wrapper function to disguise the grid points as turbines with one rotor
% point and experience almost the same calculations as the rotor points in
% the simulation.
%   Multi thread version, requires parallel computing toolbox
%   For one plot, the start-up time of the parallel toolbox can be longer
%   than the saved computational time.
Zred  = zeros(size(X));
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
    
    Zred(iGP)  = tmpM(1);
    %ZambI(iGP) = tmpM(2);
    Zeff(iGP)  = tmpM(3);
end
%U = Zeff;
Z = cat(3,Zred,Zeff);
end

function W = getWeights(X, Y, T, Dyn, TimeStep)
% GETWEIGHTS generates a matrix with weights based on the distance of the
% OPs to each other. Wind direction is taken into account, also a temporal
% decay factor is implemented. Implementation inspired by Lejeune.

distX = (X-T.States_OP(:,1)');
distY = (Y-T.States_OP(:,2)');
wPhi  = exp(- ...
    (distX.^2 + distY.^2) / (2 * Dyn.IterSigma_DW^2));
phi = (wPhi * T.States_WF(:,2))./sum(wPhi,2);


phiW = angSOWFA2world(phi);

distDW = ...
    cos(phiW).*(X-T.States_OP(:,1)') + ...
    sin(phiW).*(Y-T.States_OP(:,2)');

distCW = ...
    sin(phiW).*(X-T.States_OP(:,1)') - ...
    cos(phiW).*(Y-T.States_OP(:,2)');

W = exp(...
    - distDW.^2 / (2 * Dyn.IterSigma_DW^2) - ...
    distCW.^2 / (2 * Dyn.IterSigma_CW^2) );

W = W .* repmat(...
    exp(-((0:T.nOP-1)*TimeStep).^2 / (2 * Dyn.IterSigma_time^2)),...
    1,T.nT);
end

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

% f.Units               = 'centimeters';
% f.Position(3)         = 16.1; % A4 line width
% f.PaperPositionMode   = 'auto';
% exportgraphics(gcf,'9T_HorFlow_t300.pdf','ContentType','vector')