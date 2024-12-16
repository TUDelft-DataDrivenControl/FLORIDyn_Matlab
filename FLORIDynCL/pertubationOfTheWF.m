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

function T = pertubationOfTheWF(T,Wind)
%% pertubationOfTheWF adds noise to the entire wind field state
% Velocity
if Wind.Pertubation.Vel
    T.States_WF(:,1) = T.States_WF(:,1) + ...
        Wind.Pertubation.VelSigma * randn(T.nOP*T.nT,1);
end
% Direction
if Wind.Pertubation.Dir
    T.States_WF(:,2) = T.States_WF(:,2) + ...
        Wind.Pertubation.DirSigma * randn(T.nOP*T.nT,1);
end
% Turbulence intensity
if Wind.Pertubation.TI
    T.States_WF(:,3) = T.States_WF(:,3) + ...
        Wind.Pertubation.TISigma * randn(T.nOP*T.nT,1);
end

end