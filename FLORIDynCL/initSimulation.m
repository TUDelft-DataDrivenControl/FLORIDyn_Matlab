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

function T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS)
% Initialize the simulation or load an initialized state
switch lower(Sim.Init)
    case 'init'
        
        if Sim.SaveInitState
            save([Sim.PathToSim 'T_init.mat'],'T')
        end
    case 'load'
        % load 
        try
            load([Sim.PathToSim 'T_init.mat'],'-mat','T')
        catch
            warning(['Could not load initT.mat from ' Sim.PathToSim ...
                '\nWill proceed with initialized data.'])
        end
end
end