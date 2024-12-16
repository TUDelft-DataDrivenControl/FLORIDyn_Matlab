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

function createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,SimTime,mode)
%CREATEANDSTOREWEIGHTEDFLOWFIELD Summary of this function goes here
%   Detailed explanation goes here

%% 
fieldRes  = Vis.FlowField.Res;
fieldLims = Vis.FlowField.Lims;
xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);

[X,Y] = meshgrid(xAx,yAx);
%%
switch lower(mode)
    case 'dir'
        W_mat = getWeights(X(:), Y(:), T, Sim.Dyn.Dir, Sim.TimeStep);
        val = T.States_WF(:,2);
        V = reshape(W_mat./sum(W_mat,2) * val,size(X));
        folderName = [Sim.PathToSim 'Results' filesep 'FlowFieldRaw'...
            filesep 'Dir'];
    case 'vel'
        W_mat = getWeights(X(:), Y(:), T, Sim.Dyn.Vel, Sim.TimeStep);
        val = T.States_WF(:,1);
        V = reshape(W_mat./sum(W_mat,2) * val,size(X));
        folderName = [Sim.PathToSim 'Results' filesep 'FlowFieldRaw'...
            filesep 'Vel'];
    case 'eff'
        W_mat = getWeights(X(:), Y(:), T, Sim.Dyn.Vel, Sim.TimeStep);

        if Vis.FlowField.Plot.parallel
            Z = getMeasurementsPar2(X,Y,3,T.posNac(1,3),T,paramFLORIS,Wind,Vis);
        else
            Z = getMeasurements(X,Y,3,T.posNac(1,3),T,paramFLORIS,Wind,Vis);
        end
        
        val = T.States_WF(:,1);
        V = Z(:,:,1) .* reshape(W_mat./sum(W_mat,2) * val,size(X));
        folderName = [Sim.PathToSim 'Results' filesep 'FlowFieldRaw'...
            filesep 'Eff'];
    case 'op'
        V = T.States_OP(:,1:2);
        folderName = [Sim.PathToSim 'Results' filesep 'FlowFieldRaw'...
            filesep 'OP'];
end

if not(isfolder(folderName)); mkdir(folderName); end
writematrix(V,[folderName filesep num2str(SimTime) '.csv'])
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