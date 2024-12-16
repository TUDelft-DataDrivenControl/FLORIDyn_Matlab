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

% === Ramp event wind direction ===
varphi_0    = 225;
varphi_del  = 10;
t_0         = 31000;
t_1         = 31600;
t           = linspace(t_0,t_1);

varphi = varphi_del/2*(sin(pi*(t-t_0)/(t_1-t_0) - pi/2)+1)+varphi_0;


wind_dir = ...
    [0, varphi_0
    t', varphi';
    90000, varphi_0 + varphi_del];