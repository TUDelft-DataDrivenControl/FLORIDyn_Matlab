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

function [C_u, C_phi] = EnKF_calcC_uC_phi(truePos, trueWindDir, Dyn, Tpos, TnOP, TimeStep)
%EnKF_calcC calculates the weighted C matrix based on the state and the
%true OP position.
TnT = size(Tpos,1);
%% Get distances
distX = (Tpos(:,1) - truePos(:, 1)');
distY = (Tpos(:,2) - truePos(:, 2)');

%% Get wind direction
wPhi  = exp(- ...
    (distX.^2 + distY.^2) / (2 * Dyn.Dir.IterSigma_DW^2)); % DW & 
phi = (wPhi * trueWindDir)./sum(wPhi,2);

phiW = angSOWFA2world(phi);

%% Convert distances from x,y to down and crosswind
distDW = ...
    cos(phiW) .* distX + ...
    sin(phiW) .* distY;

distCW = ...
    sin(phiW) .* distX - ...
    cos(phiW) .* distY;

%% Calculate weights C_u
% Spacial decay
C_u = exp(...
    - distDW.^2 / (2 * Dyn.Vel.IterSigma_DW^2) - ...
    distCW.^2 / (2 * Dyn.Vel.IterSigma_CW^2) );
% Temporal decay
C_u = C_u .* repmat(...
    exp(-((0:TnOP-1)*TimeStep).^2 / (2 * Dyn.Vel.IterSigma_time^2)),...
    1,TnT);
% Norming
C_u = C_u./sum(C_u,2);
%% Calculate weights C_phi
% Spacial decay
C_phi = exp(...
    - distDW.^2 / (2 * Dyn.Dir.IterSigma_DW^2) - ...
    distCW.^2 / (2 * Dyn.Dir.IterSigma_CW^2) );
% Temporal Decay
C_phi = C_phi .* repmat(...
    exp(-((0:TnOP-1)*TimeStep).^2 / (2 * Dyn.Dir.IterSigma_time^2)),...
    1,TnT);
% Norming
C_phi = C_phi./sum(C_phi,2);
end

