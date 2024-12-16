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

function EnKF = EnKF_StoreEnStates(EnKF, M, T, iE)
%ENKF_STOREENSTATES Stores the measurements and states after the simulation
%run

%% Store updated wind field states
if EnKF.Vel.Correct
    EnKF.States.Vel(:,iE) = T.States_WF(:,1);
end
if EnKF.Dir.Correct
    EnKF.States.Dir(:,iE) = T.States_WF(:,2);
end
if EnKF.TI.Correct
    EnKF.States.TI(:,iE)  = T.States_WF(:,3);
end

%% Store OP states
EnKF.States_OP(:,EnKF.nStatesOP*(iE-1)+1:EnKF.nStatesOP*iE) = ...
    T.States_OP;

%% Store T states
EnKF.States_T(:,EnKF.nStatesT*(iE-1)+1:EnKF.nStatesT*iE) = ...
    T.States_T;

%% Store measurements
EnKF.M{iE} = [EnKF.M{iE}; M];
EnKF.Output.Pow(:,iE) = ...
    table2array(EnKF.M{iE}(end-T.nT+1:end,6));

%% Store Turbine interactions
for iT = 1:T.nT
    EnKF.Interaction{iE,iT} = ...
        [T.dep{iT}',T.intOPs{iT},T.weight{iT}];
end
if isfield(T,'C_Vel')
    EnKF.weightedInteractionVel{iE} = T.C_Vel;
end
if isfield(T,'C_Dir')
    EnKF.weightedInteractionDir{iE} = T.C_Dir;
end
end

