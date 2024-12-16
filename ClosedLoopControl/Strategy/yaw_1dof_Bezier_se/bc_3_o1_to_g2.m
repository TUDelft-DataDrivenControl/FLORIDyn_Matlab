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

function g2 = bc_3_o1_to_g2(o1, g1, rn)
%O1_TO_O2 Translation of a opimization value [0,1] into g2 (needed for bc
%and bc_d)
arguments
    o1 (:,1) double
    g1 (:,1) double
    rn (1,1) double
end

g2 = (1-o1).* ( 1/3*( -sqrt( rn^2 + 3*rn * g1) - rn + 3 * g1)) + ...
        o1 .* (-1/3*( -sqrt( rn^2 - 3*rn * g1) - rn - 3 * g1));

if sum(~isreal(g2))>0
    error('BC_3 error: Generated g2 value(s) are imaginary.')
end

end