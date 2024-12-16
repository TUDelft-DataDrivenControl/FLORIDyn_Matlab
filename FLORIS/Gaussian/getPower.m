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

function P = getPower(T,M,paramFLORIS,Con)
a   = T.States_T(T.StartI,1);
yaw = deg2rad(T.States_T(T.StartI,2));

Cp      = 4*a.*(1-a).^2;
ueff    = M(:,3);

if Con.tanhYaw
    P = 0.5*paramFLORIS.airDen*(T.D/2).^2....
        *pi.*Cp.*ueff.^3.* paramFLORIS.eta.*...
        cos(yaw).^paramFLORIS.p_p.*...
        (0.5*tanh((-yaw+deg2rad(Con.yawRangeMax))*50)+.5) .* ...
        (-0.5*tanh((-yaw+deg2rad(Con.yawRangeMin))*50)+.5);
else
    P = 0.5*paramFLORIS.airDen*(T.D/2).^2....
        *pi.*Cp.*ueff.^3.* paramFLORIS.eta.*...
        cos(yaw).^paramFLORIS.p_p;
end
end