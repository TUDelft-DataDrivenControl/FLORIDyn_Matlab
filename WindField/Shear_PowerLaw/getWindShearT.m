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

function shear = getWindShearT(WindShear,z_norm)
%GETWINDSHEART Return the shear factor u_eff = shear * u_referenceHight
% POWER LAW implementation
%   expects a WindShearPowerLaw.csv with the shear coefficient
% 
% ======================================================================= %
% WindShear = Holds shear coefficient and reference height
%          .z0      = reference height
%          .alpha   = shear coefficient
% z         = height(s)
% ======================================================================= %
shear = (z_norm).^WindShear.alpha;
end

