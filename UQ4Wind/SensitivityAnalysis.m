close all
clear all
clc
clearvars
rng default

global n m sim_p;  % # of turbines % Still needed?
n = 3;
m = 600;   % choice of PCEOpts.ExpDesign.NSamples for creating surrogate model
sim_p = 0;
%% Reset the Matlab Path and load essential paths & the simulation path
% Link folder with simulation data
pathToSimulation = './Simulations/2021_9T_Data/';%2021_9T_Data/';
addPaths;
addpath(genpath('./Simulations/2021_9T_Data/'));%2021_9T_Data/'))
addpath('./UQ4Wind')
addpath(genpath('C:\Users\marcusbecker\surfdrive\PhD_Surf\01_Research\07_UQLab\UQLabCore_Rel1.4.0'))
%% Run UQLab
dispstat('','init') % Can be deleted
timeEst();
uqlab

%% Forward model

Model.mHandle = @genrunFLORIDyn;
myModel = uq_createModel(Model);
Model.isVectorized = true;

numofturbines = num2str(n);
numofsamples = num2str(m);
surrogatemodelname = strcat(numofturbines,'turbines',numofsamples,'samplesPCEsurrogatemodel4sa.mat');
%% Prior distribution of model parameters 

Input.Marginals(1).Name = 'alpha'; 
Input.Marginals(1).Type = 'Uniform'; 
%Input.Marginals(1).Parameters = [.5 3.5];%[2 2.5];
Input.Marginals(1).Parameters = [0.8704    1.3056];


Input.Marginals(2).Name = 'beta';
Input.Marginals(2).Type = 'Uniform'; 
%Input.Marginals(2).Parameters = [0.025 0.25];%[0.125 0.175];
Input.Marginals(2).Parameters = [0.1776    0.2664];


Input.Marginals(3).Name = 'ka';
Input.Marginals(3).Type = 'Uniform'; 
%Input.Marginals(3).Parameters = [0.1,0.6];%[0.35 0.4];
Input.Marginals(3).Parameters = [0.4296    0.6444];


Input.Marginals(4).Name = 'kb'; 
Input.Marginals(4).Type = 'Uniform'; 
%Input.Marginals(4).Parameters = [0.001 0.01];%[0.003 0.004];
Input.Marginals(4).Parameters = [-0.0010 -0.0007];

Input.Marginals(5).Name = 'obja';
Input.Marginals(5).Type = 'Uniform'; 
%Input.Marginals(5).Parameters = [0.7 0.75];
Input.Marginals(5).Parameters = [6.2720    9.4080];


Input.Marginals(6).Name = 'objb';
Input.Marginals(6).Type = 'Uniform'; 
%Input.Marginals(6).Parameters = [0.8 0.85];
Input.Marginals(6).Parameters = [3.6560    5.4840];


Input.Marginals(7).Name = 'objc'; 
Input.Marginals(7).Type = 'Uniform'; 
%Input.Marginals(7).Parameters = [0.03 0.035];
Input.Marginals(7).Parameters = [0.344 0.516];


Input.Marginals(8).Name = 'objd';
Input.Marginals(8).Type = 'Uniform'; 
%Input.Marginals(8).Parameters = [-0.6 -0.03];%[-0.45 -0.15];
Input.Marginals(8).Parameters = [-0.2952 -0.1968];


Input.Marginals(9).Name = 'd'; 
Input.Marginals(9).Type = 'Uniform';
Input.Marginals(9).Parameters = [0.8 1.2]; 
 

Input.Marginals(10).Name = 'eta';
Input.Marginals(10).Type = 'Uniform';
Input.Marginals(10).Parameters = [0.864 1.296];      

Input.Marginals(11).Name = 'pp'; 
Input.Marginals(11).Type = 'Uniform';  
Input.Marginals(11).Parameters = [1.76 2.64];       % pp = 1.50

myInput = uq_createInput(Input) ;
uq_print(myInput);
uq_display(myInput);

%% Surrogate model

% Surrogate model options

% Switch for sensitivity analysis with the full FLORIDyn model or with the surrogate model
Sensitivity_full = 0; % 0: use surrogate model (PCE); 1: run full model (Computationally expensive!)

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

SobolOpts.Sobol.Order = 2;

%SobolOpts.Sobol.SampleSize = 1e4;

mySobolAnalysis = uq_createAnalysis(SobolOpts);

%% Post-processing
plotTurbineSensitivity;
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