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

function [turbineProperties] = estimator_dtu10mw()
% NEW (information straight from SOWFA):
load('cpInterp_DTU10MW_FAST.mat')

% Define turbine properties
turbineProperties = struct(...
    'GearboxRatio',50.0,... % Gearbox ratio [-]
    'InertiaTotal',1.409969209E+08,... % Total inertia
    'RotorRadius',89.2,... % Rotor radius [m]
    'CpFun',cpInterpolant,... % Cp interpolant for TSR and blade pitch
    'FluidDensity',1.23,...
    'GearboxEff',1); % Fluid density
end