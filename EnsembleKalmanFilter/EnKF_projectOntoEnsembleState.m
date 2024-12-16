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

function EnKF = EnKF_projectOntoEnsembleState(EnKF,stateIndex)

for iE = 1:EnKF.nE
    switch stateIndex
        case 1
            EnKF.States.Vel(:,iE) = ...
                linsolve(EnKF.Vel.W{iE},EnKF.States.Vel(:,iE));
        case 2
            EnKF.States.Dir(:,iE) = EnKF.Dir.W{iE}\EnKF.States.Dir(:,iE);
        case 3
            error('Not implemented yet')
    end
end
end