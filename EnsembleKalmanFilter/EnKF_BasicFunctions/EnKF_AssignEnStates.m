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

function T = EnKF_AssignEnStates(EnKF,T,iE)
%ENKF_ASSIGNENSTATES Assign States if the wind field state is being
%corrected for the next simulation run.

%% Ensemble States
if EnKF.Vel.Correct
    T.States_WF(:,1) = EnKF.States.Vel(:,iE);
end
if EnKF.Dir.Correct
    T.States_WF(:,2) = EnKF.States.Dir(:,iE);
end
if EnKF.TI.Correct
    T.States_WF(:,3) = EnKF.States.TI(:,iE);
end

%% OP States
T.States_OP = ...
    EnKF.States_OP(:,EnKF.nStatesOP*(iE-1)+1:EnKF.nStatesOP*iE);

%% T States
T.States_T = ...
    EnKF.States_T(:,EnKF.nStatesT*(iE-1)+1:EnKF.nStatesT*iE);
end

% === TODO ===
% Uncorrected states are just passed on now, from one ensemble to the next.
% This mostly concerns the uncorrected WF states (e.g. Ambient TI)
%   No issue if the value does not change anyways.