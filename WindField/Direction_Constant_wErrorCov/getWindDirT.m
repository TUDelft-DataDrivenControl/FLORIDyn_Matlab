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

function phi = getWindDirT(WindDir,iT,~)
%GETWINDDIRT Return wind direction in SOWFA-deg for the requested
%turbine(s)
%
% ======= Input ======
% WindDir.Data   = float, wind direction
% WindDir.ColSig = nT x nT, col(Covariance Matrix)
% iT        = Index/Indeces of the turbines

phi = ones(size(iT))*WindDir.Data;
phi = phi + (randn(1,length(phi))*WindDir.CholSig)';
end

