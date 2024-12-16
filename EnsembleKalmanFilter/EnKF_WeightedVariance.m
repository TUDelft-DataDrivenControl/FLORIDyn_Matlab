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

function [trueOPpos, meanVal, varVal, stdVal] = EnKF_WeightedVariance(EnKF, T, Vis, Sim, iState)
%ENKF_WEIGHTEDVARIANCE calculates the mean and variance of multiple
%ensembles based on a spacial-temporal weighting of the value.
%   In a first step the mean loaction of the 'true' OPs is determined. Then
%   the states of each ensemble are used to find out what the weighted
%   value at the location of the true state is. For that, first the weights
%   need to be determined by distance and OP age, then the value can be
%   calculated. Finally, the value for the state is determined in each
%   ensemble and a mean and variance can be calculated across all
%   ensembles.
%% Step 0: Load appropriate data
switch iState
    case 1
        Dyn = Sim.Dyn.Vel;
        Val = EnKF.States.Vel;
    case 2
        Dyn = Sim.Dyn.Dir;
        Val = EnKF.States.Dir;
    otherwise
        error('Only iState==1 or 2 defined. 1 -> Velocity, 2 -> Direction')
end
valTrue = zeros(size(Val));
%% Step 1: Calculate location of 'true' OPs
trueOPpos = [...
    mean(EnKF.States_OP(:,1:length(T.Names_OP):end),2),...
    mean(EnKF.States_OP(:,2:length(T.Names_OP):end),2)];

%% Step 2: Cycle through all ensembles and ...
for iE = 1:EnKF.nE
    T = EnKF_AssignEnStates(EnKF,T,iE);

%% Step 2.1: Calculate the weights based on distance and OP age
    W = getWeights(trueOPpos(:,1), trueOPpos(:,2), T, Dyn, Sim.TimeStep);
    W = W./sum(W,2);
    
%% Step 2.2: Calculate the values at the true OPs based on the weights
    valTrue(:,iE) = W * Val(:,iE);

end
%% Step 3: Collect all values across all ensembles for the true OPs
meanVal = mean(valTrue,2);
varVal = var(valTrue,[],2);
stdVal = std(valTrue,[],2);
end


%% Additional functions
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