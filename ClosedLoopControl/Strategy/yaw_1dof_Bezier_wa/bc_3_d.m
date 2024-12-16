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

function bc_d = bc_3_d(t_n,g0,g1,g2,g3)
%BC_3_D is the derivative of a cubic bezier-curve, with t_n \in [0,1].
%   Returns a matrix of the size (timesteps x turbines)
arguments
    t_n (:,1) double
    g0 (1,:) double
    g1 (1,:) double
    g2 (1,:) double
    g3 (1,:) double
end

bc_d = 3*(1-t_n).^2  * (g1 - g0)  +...
    6*(1-t_n).* t_n * (g2 - g1) +...
    3*t_n.^2 * (g3 - g2);

end

