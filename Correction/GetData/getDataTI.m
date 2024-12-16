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

function TI = getDataTI(Wind,T,SimTime)
%GETDATATI retrieves the data for the ambient turbulence intensity

switch Wind.Input_TI
    case 'EnKF_InterpTurbine'
        TI = getWindTiT_EnKF(Wind.TI,(1:T.nT)',SimTime);
    case 'EnKF_ZOH'
        TI = T.States_WF(T.StartI,3);
    case 'EnKF_RW'
        TI = T.States_WF(T.StartI,3);
        TI = TI + (randn(1,length(phi))*Wind.TI.CholSig)';
    case 'CLC_weighted_ZOH'
        % Apply weighted states on the turbine
        TI = T.C_TI*T.States_WF(:,3);
    otherwise
        TI = getWindTiT(Wind.TI,(1:T.nT)',SimTime);
end
end

