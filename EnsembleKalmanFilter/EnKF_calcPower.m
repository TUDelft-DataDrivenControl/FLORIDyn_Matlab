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

function [P,EnKF] = EnKF_calcPower(C_u, EnKF, TD, paramFLORIS, TStartI)
%ENKF_CALCPOWER Calculates the power of the turbines
%   This is done as the wind field states have been projected on a common
%   OP location, at which we correct the states.

% Free wind speed from the ensembles
Vel = EnKF.States.Vel;
nT  = length(TD);
red = zeros(nT, EnKF.nE);
yaw = deg2rad(EnKF.States_T(TStartI,2:EnKF.nStatesT:end));

for iE = 1:EnKF.nE
    red(:,iE) = EnKF.M{iE}.("Foreign Reduction [%]")(end-nT+1:end)...
        * 10^(-2);
end

ueff = red .* (C_u * Vel);
a    = 1/3;
Cp   = 4*a.*(1-a).^2;
P    = 0.5* paramFLORIS.airDen*(TD/2).^2....
        .* pi .* Cp .* ueff.^3 .* paramFLORIS.eta * 10^(-6) .* ...
        cos(yaw).^paramFLORIS.p_p ;

% Relace the measurements for the new ones
EnKF = EnKF_overwritePowerM(EnKF,nT,P);

end

