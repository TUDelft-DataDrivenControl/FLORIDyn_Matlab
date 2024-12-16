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


close all
clear all
clc
clearvars
rng default

global n m sim_p;  % # of turbines % Still needed?
n = 9;
m = 400;   % choice of PCEOpts.ExpDesign.NSamples for creating surrogate model
sim_p = 0;
%% Reset the Matlab Path and load essential paths & the simulation path
% Link folder with simulation data
pathToSimulation = './Simulations/2021_9T_Data/';%2021_9T_Data/';
addPaths;
addpath(genpath('./Simulations/2021_9T_Data/'));%2021_9T_Data/'))
addpath('./UQ_EnKF')
addpath(genpath('C:\Users\marcusbecker\Downloads\UQLab_Rel2.0.0\UQLab_Rel2.0.0'))

%% Run UQLab
%dispstat('','init') % Can be deleted
%timeEst();
uqlab

%% Forward model

Model.mHandle = @genrunFLORIDynEnKF;
myModel = uq_createModel(Model);
Model.isVectorized = true;

numofturbines = num2str(n);
numofsamples = num2str(m);
path2Storage = 'O:\TUDelft\EnKF_Sensitivity_Data\'; %'/Volumes/MarcusSSD/TUDelft/EnKF_Sensitivity_Data/';
surrogatemodelname = strcat(path2Storage,numofturbines,'turbines',numofsamples,'samplesPCEsurrogatemodel4sa.mat');
%% Prior distribution of model parameters 

Input.Marginals(1).Name = 'velWeight_{dw}'; 
Input.Marginals(1).Type = 'Uniform'; 
Input.Marginals(1).Parameters = [50    600];


Input.Marginals(2).Name = 'velWeight_{cw}';
Input.Marginals(2).Type = 'Uniform'; 
Input.Marginals(2).Parameters = [25    250];


Input.Marginals(3).Name = 'velWeight_t';
Input.Marginals(3).Type = 'Uniform'; 
Input.Marginals(3).Parameters = [50    500];


Input.Marginals(4).Name = 'dirWeight_{dw}'; 
Input.Marginals(4).Type = 'Uniform'; 
Input.Marginals(4).Parameters = [200    1000];

Input.Marginals(5).Name = 'dirWeight_{cw}';
Input.Marginals(5).Type = 'Uniform'; 
Input.Marginals(5).Parameters = [200    1000];


Input.Marginals(6).Name = 'dirWeight_t';
Input.Marginals(6).Type = 'Uniform'; 
Input.Marginals(6).Parameters = [25    150];


Input.Marginals(7).Name = 'n_E'; 
Input.Marginals(7).Type = 'Uniform';
Input.Marginals(7).Parameters = [25    150];

myInput = uq_createInput(Input) ;
uq_print(myInput);
uq_display(myInput);

%% Surrogate model

% Surrogate model options

% Switch for sensitivity analysis with the full FLORIDyn model or with the surrogate model
Sensitivity_full = 1; % 0: use surrogate model (PCE); 1: run full model (Computationally expensive!)

% If Sensitivity_full = 0, we need to specify options for creating or loading a surrogate model
Surrogate_model = 1; % 0: use a stored PCE surrogate model, 1: create surrogate model


if (Sensitivity_full == 0) % create or load a PCE surrogate model to be used
    if (Surrogate_model == 0)
        disp('using stored surrogate model');
        load(surrogatemodelname) % saved surrogates as .mat files
        SobolOpts.Model = myPCE;
    elseif (Surrogate_model == 1)
        disp('creating surrogate model');
        % use prior also as input uncertainties
        PCEOpts.Type = 'Metamodel';
        PCEOpts.MetaType = 'PCE';
        PCEOpts.FullModel = myModel;
        PCEOpts.Degree = 1:4;
        PCEOpts.TruncOptions.qNorm = 0.7;
        PCEOpts.TruncOptions.MaxInteraction = 2;                                     
        PCEOpts.ExpDesign.NSamples = m;
        PCEOpts.ExpDesign.Sampling = 'LHS';
        myPCE = uq_createModel(PCEOpts);
        
        %str1= 'D:\FLORIDyn3\FLORIDyn_Matlab\surrogatemodel\SensitivityAnalysis\';
        %numofturbines = num2str(n);
        %numofsamples = num2str(m);
        %str2 = strcat(numofturbines,'turbines',numofsamples,'PCEsurrogatemodel');
        %surrogatemodelname = strcat(str1,str2);
        save(surrogatemodelname,'myPCE');
          
        SobolOpts.Model = myPCE;
    end
else % do SA with full model
    SobolOpts.Model = myModel;
end

%% Sensitivity analysis

SobolOpts.Type        = 'Sensitivity';
SobolOpts.Method      = 'Sobol';

SobolOpts.Sobol.Order = 1;

%SobolOpts.Sobol.SampleSize = 1e4;

mySobolAnalysis = uq_createAnalysis(SobolOpts);

%% Post-processing
%plotTurbineSensitivity;
% figure
% bar(mySobolAnalysis.Results.Total)
% 
% %%
% figure()
% a = area(mySobolAnalysis.Results.Total','LineWidth',1.5,'LineStyle','-');
% newcolor = [[255 0 0]/255;...
%             [255 135 0]/255;...
%             [255 211 0]/255;...
%             [222 255 10]/255;...
%             [161 255 10]/255;...
%             [10 255 153]/255;...
%             [10 239 255]/255;...
%             [20 125 245]/255;...
%             [2 62 125]/255;...
%             [141 153 174]/255;...  
%             [45 35 46]/255   ];
%         colororder(newcolor)   
% ylim([0 1])
% legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','\theta_8','\theta_9','\theta_{10}','\theta_{11}','Orientation','horizontal','Location','southoutside','Box','off','FontSize', 14)
% xlim([502 753])
% grid on
% xticklabels({200 400 600 800 1000})
% %ylabel('time [s]')
% xlabel('time [s]')
% ylabel('S^{Total}_i')
disp('Done.')