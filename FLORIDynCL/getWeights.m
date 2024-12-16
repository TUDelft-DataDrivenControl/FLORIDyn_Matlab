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

function W = getWeights(X, Y, T, Dyn, TimeStep)
% GETWEIGHTS generates a matrix with weights based on the distance of the
% OPs to each other. Wind direction is taken into account, also a temporal
% decay factor is implemented. Implementation inspired by Lejeune.

max_sig = max(Dyn.IterSigma_DW,Dyn.IterSigma_CW);          % Addition to reduce calculations

distX = (X-T.States_OP(:,1)');
distY = (Y-T.States_OP(:,2)');

consider = or(abs(distX)<3*max_sig, abs(distY)<3*max_sig); % Addition to reduce calculations

wPhi = zeros(distX);

wPhi(consider)  = exp(- ...
    (distX(consider).^2 + distY(consider).^2) / (2 * Dyn.IterSigma_DW^2));

phi = (wPhi * T.States_WF(:,2))./sum(wPhi,2);
                                    
phiW = angSOWFA2world(phi);

distDW = ...
    cos(phiW).*distX(consider) + ...
    sin(phiW).*distY(consider);

distCW = ...
    sin(phiW).*distX(consider) - ...
    cos(phiW).*distY(consider);

W = zeros(size(distX));

W(consider) = exp(...
    - distDW.^2 / (2 * Dyn.IterSigma_DW^2) - ...
    distCW.^2 / (2 * Dyn.IterSigma_CW^2) );

W = W .* repmat(...
    exp(-((0:T.nOP-1)*TimeStep).^2 / (2 * Dyn.IterSigma_time^2)),...
    1,T.nT);
end
