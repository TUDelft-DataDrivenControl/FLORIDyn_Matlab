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

function tr = bc_3_tr(t_n,g0,g1,g2,g3)
%BC_3 is the trajectory of a cubic bezier-curve, with t_n \in [0,1].
%   Returns a matrix of the size (timesteps x turbines)
arguments
    t_n (:,1) double
    g0 (1,:) double
    g1 (1,:) double
    g2 (1,:) double
    g3 (1,:) double
end

tr = (1-t_n).^3*g0 + 3*(1-t_n).^2 .* t_n * g1 + ...
    3*(1-t_n).* t_n.^2 * g2 + t_n.^3 * g3;

end

