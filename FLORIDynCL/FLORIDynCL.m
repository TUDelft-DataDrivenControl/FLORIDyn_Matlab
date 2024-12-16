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

function [T,Mt,Vis,Mint] = FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS)
%FLORIDYNCL Summary of this function goes here
%   
% ======================================================================= %
% INPUTS
% T     : Turbine data and OP states
% Wind  : Data for the environment (wind speed, direction, ...)
% Sim   : Simulation parameters and settings
% Con   : Settings regarding the control design
% Vis   : Settings regarding the Visualization
% paramFLORIDyn : Parameter for FLORIDyn
% paramFLORIS   : Parameter for FLORIS
% ======================================================================= %
% OUTPUTS
% T := Simulation state (OP states, Turbine states, wind field states(OPs))
% M := Measurements from the simulation (Power, tbd)
% Vis := Visualization data
% Mint := Interaction matrices for the turbines
% ======================================================================= %
%% Run simulation
% T_red, T_addedTI, T_Ueff
M = zeros(Sim.nSimSteps*T.nT,6);
M(:,1) = 1;
M_int = cell(Sim.nSimSteps,1);

Vis.StoreFlowField.Vel      = false;
Vis.StoreFlowField.Dir      = false;
Vis.StoreFlowField.Eff      = false;
Vis.StoreFlowField.OP       = false;
Vis.FlowField.Plot.parallel = true;

SimTime = Sim.StartTime;
for it = 1:Sim.nSimSteps
    Sim.SimStep = it;
    % ========== PREDICTION ==========
    % Iterate OPs and states
    T = iterateOPs(T,Sim,paramFLORIS,paramFLORIDyn);
    
    % ========== Wind Field Pertubation
    T = pertubationOfTheWF(T,Wind);
    
    % ========== Get FLORIS reductions
    %   ======== Find groups
    T.dep = findTurbineGroups(T, paramFLORIDyn);
    
    %   ======= Interpolate OPs
    T.intOPs = interpolateOPs(T);
    
    %   ======= Set up temp Wind farms & run FLORIS
    [tmpM,T] = setUpTmpWFAndRun(T,paramFLORIS,Wind);
    M((it-1)*T.nT+1:it*T.nT,2:4)    = tmpM;
    M((it-1)*T.nT+1:it*T.nT,1)      = SimTime;
    T.States_T(T.StartI,3)          = tmpM(:,2);
    M_int{it}                       = T.red_arr;
    
    % ========== Get Wind field variables & correct values
    [T,Wind] = correctVel(T,Wind,SimTime,paramFLORIS,tmpM);
    T = correctDir(T,Wind,SimTime);
    T = correctTi(T,Wind,SimTime);
    
    % Save free wind speed as measurement
    M((it-1)*T.nT+1:it*T.nT,5) = T.States_WF(T.StartI,1);
    
    % ========== Get Control settings
    % Yaw is read in orientation angles
    %T.States_T(T.StartI,1) = getAxInd(Con.YawData,(1:T.nT)',SimTime);     % PLACEHOLDER for axial induction control 
    T.States_T(T.StartI,2) = (T.States_WF(T.StartI,2) - ...
        getYaw(Con.YawData,(1:T.nT)',SimTime)');
    
    % ========== Calculate Power
    P = getPower(T,tmpM,paramFLORIS,Con);
    
    M((it-1)*T.nT+1:it*T.nT,6) = P;
    
    % ========== Live Plotting
    if Vis.FlowField.Plot.Online
        Vis = PlotFlowField(T,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS,...
            SimTime);
    end
    if Vis.FlowField.Error.Online
        if sum(Vis.FlowField.Error.Steps==SimTime)==1
            errorVTK([Vis.FlowField.Error.ValidationPath ...
                num2str(SimTime) '\' Vis.FlowField.Error.SnapName],...
                [num2str(SimTime) 'FlowError'],...
                T,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS)
        end
    end
    if Vis.StoreFlowField.Vel
        createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
            SimTime,'Vel')
    end
    if Vis.StoreFlowField.Dir
        createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
            SimTime,'Dir')
    end
    if Vis.StoreFlowField.Eff
        createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
            SimTime,'Eff')
    end
    if Vis.StoreFlowField.OP
        createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
            SimTime,'OP')
    end
    % ========== Increase time
    SimTime = SimTime + Sim.TimeStep;
    
    % ==== Store the centerline data ====
    if Vis.CL.Store
        storeCLInfo(T,Sim,it);
    end
end

if Vis.FlowField.Plot.Online
    Vis.Film.InProgress = false;
    % Write saved frames as a movie to the simulation folder
    v = VideoWriter(Vis.Film.MovFileEffU);
    v.FrameRate = 12;
    v.Quality   = 100;
    open(v)
    writeVideo(v,Vis.Film.FrmFileEffU)
    close(v)
end

% Convert M to table and scale measurements accordingly
Mt = array2table(M*diag([1;100;100;1;1;10^(-6)]),...
    'VariableNames',{'Time [s]','Foreign Reduction [%]','Added turbulence [%]',...
    'Effective Wind Speed [ms^-1]','Free Wind Speed [ms^-1]', ...
    'Power generated [MW]'});
Mint = [Mt.("Time [s]") cell2mat(M_int)];
end
