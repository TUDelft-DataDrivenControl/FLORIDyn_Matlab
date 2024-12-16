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

%% Run UQLab

uqlab

%% Forward model
Model.mHandle = @genrunFLORIDyn;
myModel = uq_createModel(Model);
Model.isVectorized = true;

global n;                                             % # of turbines
m = 200;                                            % choice of PCEOpts.ExpDesign.NSamples for creating surrogate model

numofturbines = num2str(n);
numofsamples = num2str(m);
%surrogatemodelname = strcat(numofturbines,'turbines',numofsamples,'samplesPCEsurrogatemodel4bc.mat');
surrogatemodelname = 'C:\Users\kriek\Downloads\Group_presentation/3t_PCE.mat';
%% Experimental data

exp = load('C:\Users\kriek\Downloads\Group_presentation/3tSOWFA.mat');
exp_f = [exp.p(:,1);exp.p(:,2);exp.p(:,3)]';   %SOWFA data


myData.y = exp_f; % (m)
myData.Name = 'Power';

%% Prior distribution of model parameters 

Input.Marginals(1).Name = 'alpha'; 
Input.Marginals(1).Type = 'Uniform'; 
Input.Marginals(1).Parameters = [2 2.5];


Input.Marginals(2).Name = 'beta';
Input.Marginals(2).Type = 'Uniform'; 
Input.Marginals(2).Parameters = [0.125 0.175];


Input.Marginals(3).Name = 'ka';
Input.Marginals(3).Type = 'Uniform'; 
Input.Marginals(3).Parameters = [0.35 0.4];


Input.Marginals(4).Name = 'kb'; 
Input.Marginals(4).Type = 'Uniform'; 
Input.Marginals(4).Parameters = [0.003 0.004];

Input.Marginals(5).Name = 'obja';
Input.Marginals(5).Type = 'Uniform'; 
Input.Marginals(5).Parameters = [0.7 0.75];


Input.Marginals(6).Name = 'objb';
Input.Marginals(6).Type = 'Uniform'; 
Input.Marginals(6).Parameters = [0.8 0.85];


Input.Marginals(7).Name = 'objc'; 
Input.Marginals(7).Type = 'Uniform'; 
Input.Marginals(7).Parameters = [0.03 0.035];


Input.Marginals(8).Name = 'objd';
Input.Marginals(8).Type = 'Uniform'; 
Input.Marginals(8).Parameters = [-0.45 -0.15];


Input.Marginals(9).Name = 'd'; 
Input.Marginals(9).Type = 'Uniform';
Input.Marginals(9).Parameters = [0.7 1.3]; 
 

Input.Marginals(10).Name = 'eta';
Input.Marginals(10).Type = 'Uniform';
Input.Marginals(10).Parameters = [1 1.2];      

Input.Marginals(11).Name = 'pp'; 
Input.Marginals(11).Type = 'Uniform';  
Input.Marginals(11).Parameters = [1.3 1.6];       % pp = 1.50

myInput = uq_createInput(Input) ;
uq_print(myInput);
uq_display(myInput);

%% Discrepancy

DiscrepancyOptsKnown.Type = 'Gaussian'; % known
DiscrepancyOptsKnown.Parameters = 1e-2;

%% Surrogate model
% Surrogate model options

% Switch for Bayesian analysis with the full FLORIDyn model or with the surrogate model
Bayesian_full = 0; % 0: use surrogate model (PCE); 1: run full model (Computationally expensive!)

% If Bayesian_full = 0, we need to specify options for creating or loading a surrogate model
Surrogate_model = 0; % 0: Uses a stored PCE surrogate model, 1: create surrogate model

if (Bayesian_full == 0) % create or load a PCE surrogate model to be used
    if (Surrogate_model == 0)
        disp('using stored surrogate model');
        load(surrogatemodelname) % saved surrogates as .mat files
        BayesOpts.ForwardModel.Model = mySurrogateModel;
    elseif (Surrogate_model == 1)
        disp('creating surrogate model');
        % use prior also as input uncertainties
        MetaOpts.Type = 'Metamodel';
        MetaOpts.MetaType = 'PCE';
        MetaOpts.ExpDesign.NSamples = m;
        mySurrogateModel = uq_createModel(MetaOpts);
        
        save(surrogatemodelname,'mySurrogateModel');
        
        BayesOpts.ForwardModel.Model = mySurrogateModel;
    end
else % do Bayesian calibration with full model
    BayesOpts.ForwardModel.Model = myModel;
end

%% Solver options

Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES'; % MH 
Solver.MCMC.Steps = 150;
Solver.MCMC.NChains = 1000;

%% Bayes options

BayesOpts.Type = 'Inversion';
BayesOpts.Data = myData;
BayesOpts.Prior = myInput;
BayesOpts.Solver = Solver;
BayesOpts.Discrepancy = DiscrepancyOptsKnown;

myBayesianAnalysis = uq_createAnalysis(BayesOpts);

%% Print analysis

uq_print(myBayesianAnalysis)

%% Post-processing

uq_display(myBayesianAnalysis)

uq_display(myBayesianAnalysis, 'trace', 'all')

%% Cross-validation

%Sim.UQ = false;

