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

function phi = getDataDir(Wind,T,SimTime)
%GETDATADIR Reads wind data and returns the current phi for all turbines

switch Wind.Input_Dir
    case 'RW_with_Mean'
        phi = getWindDirT(T.States_WF(T.StartI,2),Wind.Dir);
    case 'EnKF_ZOH'
        phi = T.States_WF(T.StartI,2);
    case 'EnKF_RW'
        phi = T.States_WF(T.StartI,2);
        phi = phi + (randn(1,length(phi))*Wind.Dir.CholSig)';
    case 'EnKF_InterpTurbine'
        phi = getWindDirT_EnKF(Wind.Dir,(1:T.nT)',SimTime);
    case 'CLC_weighted_ZOH'
        % Apply weighted states on the turbine
        phi = T.C_Dir*T.States_WF(:,2);
    otherwise
        phi = getWindDirT(Wind.Dir,(1:T.nT)',SimTime);
end

end

