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

function T = iterateOPs(T,Sim,paramFLORIS,paramFLORIDyn)
%ITERATEOPS propagates the OPs downstream and deletes the oldest ones.


%% Save turbine OPs
tmpOPStates = T.States_OP(T.StartI,:);
tmpTStates  = T.States_T(T.StartI,:);
tmpWFSTates = T.States_WF(T.StartI,:);
%% Shift states
% Calculate downwind step and apply to wake coordinate system
step_dw = Sim.TimeStep * T.States_WF(:,1) * Sim.Dyn.Advection;
T.States_OP(:,4) = T.States_OP(:,4) + step_dw;

% Calculate crosswind step and apply to wake coordinate system
%   TODO If turbines with different diameters are used, this has to be run
%   individually (in parallel) for all turbines.
deflection  = Centerline(...
    T.States_OP,T.States_T,T.States_WF,paramFLORIS,T.D(1));
step_cw     = (deflection - T.States_OP(:,5:6))*1;

T.States_OP(:,5:6) = deflection*1;

% Apply dw and cw step to the world coordinate system
phiW = angSOWFA2world(T.States_WF(:,2));

T.States_OP(:,1) = T.States_OP(:,1) + ...
    cos(phiW) .* step_dw - ...
    sin(phiW) .* step_cw(:,1);
T.States_OP(:,2) = T.States_OP(:,2) + ...
    sin(phiW) .* step_dw + ...
    cos(phiW) .* step_cw(:,1);
T.States_OP(:,3) = T.States_OP(:,3) + step_cw(:,2);

%% Circshift & Init first OPs
%   OPs
T.States_OP = circshift(T.States_OP,1);
T.States_OP(T.StartI,:) = tmpOPStates;
%   Turbines
T.States_T  = circshift(T.States_T,1);
T.States_T(T.StartI,:)  = tmpTStates;
%   Wind Field
% Getting weights for averaging
wNew = Sim.Dyn.OPiterWeights(1);
wOld = Sim.Dyn.OPiterWeights(2);
T.States_WF = wOld*T.States_WF + wNew*circshift(T.States_WF,1);
T.States_WF(T.StartI,:) = tmpWFSTates;

%% Check if OPs are in order
for iT=1:T.nT
    [~,indOP] = sort(T.States_OP(T.StartI(iT)+(0:(T.nOP-1)),4));
    if ~issorted(indOP)
        warning(['OPs overtaking, consider increasing the weight on ' ...
            'the old wind field state in setup > OP / wind field propagation.'])
        T.States_OP(T.StartI(iT)+(0:(T.nOP-1)),:) = ...
            T.States_OP(T.StartI(iT) + indOP - 1,:);
        T.States_T(T.StartI(iT)+(0:(T.nOP-1)),:) = ...
            T.States_T(T.StartI(iT) + indOP - 1,:);
        T.States_WF(T.StartI(iT)+(0:(T.nOP-1)),:) = ...
            T.States_WF(T.StartI(iT) + indOP - 1,:);
    end
end

end

