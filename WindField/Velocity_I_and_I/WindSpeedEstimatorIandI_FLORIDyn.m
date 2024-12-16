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

function [V_out,WSE] = WindSpeedEstimatorIandI_FLORIDyn(WSE, Rotor_Speed, Blade_pitch, Gen_Torque,yaw,p_p)
%% Extended I&I Wind speed estimator
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Reference: 
% Liu Y., Pamososuryo A., Ferrari M.G. R., van Wingerden J.W., 
% The Immersion and Invariance Wind Speed Estimator Revisited and New Results, 
% IEEE Control Systems Letters, 2021.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% The script is based on the work of Yichao Liu and Jean Gonzales Silva.
% Modified by Marcus Becker in 22/07/2021.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%% Inputs
% WSE       : Wind Speed Estimator States and Parameters
%    .Ee    : State - Integrated error between est. rot. speed and measured
%    .V     : State - Estimated wind speed [m s^-1]
%    .omega : State - Estimated rotor speed [rad s^-1]
%    .beta  : Integral gain [-]
%    .gamma : Proportional gain (>0) [-]
%    .T_prop: Turbine properties
%           .gearboxratio   : Gearbox ratio [-]
%           .inertiaTotal   : Total inertia 
%           .rotorRadius    : Rotor radius [m]
%           .cpFun          : Cp interpolant for TSR and blade pitch
%           .fluidDensity   : Fluid density
% Rotor_Speed   : Rotor speed [rpm]
% Blade_pitch   : Blade pitch angle [deg]
% Gen_Torque    : Generator torque [Nm]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%% Import variables from WSE
fluidDensity    = WSE.T_prop.FluidDensity;  % Fluid density [kg/m3]
rotorRadius     = WSE.T_prop.RotorRadius;   % Rotor radius [m]
rotorArea       = pi*rotorRadius^2;         % Rotor swept surface area [m2]
gamma           = WSE.gamma;                % Estimator gain
beta            = WSE.beta;                 % Estimator gain
gbRatio         = WSE.T_prop.GearboxRatio;  % Gearbox ratio
inertTot        = WSE.T_prop.InertiaTotal;  % Inertia
dt              = WSE.dt_SOWF;              % Timestep_SOWFA
GBEfficiency    = WSE.T_prop.GearboxEff;    % Gearbox efficiency
%% Calculate aerodynamic torque
%   Estimated tip speed ratio [-]
tipSpeedRatio   = (Rotor_Speed * pi * rotorRadius)./(WSE.V * 30); 
%   Power coefficient [-]
Cp = max(WSE.T_prop.CpFun(tipSpeedRatio,Blade_pitch),0)*fluidDensity; 
%   Correct for yaw angle
Cp = Cp.*cos(yaw).^p_p;
%% Estimate wind speed
if isnan(Cp)
    disp(['Cp is out of the region of operation: TSR=' ...
        num2str(tipSpeedRatio) ', Pitch=' num2str(Blade_pitch) ' deg.'])
    disp('Assuming windSpeed to be equal to the past time instant.')
else
    aerodynamicTorque = 0.5 * fluidDensity * rotorArea *...
        ((WSE.V.^3)./(Rotor_Speed* pi/30)) .* Cp; % Torque [Nm]
    
    % Saturate torque to non-negative numbers
    aerodynamicTorque = max(aerodynamicTorque, 0.0); 

    % Update estimator state and wind speed estimate (YICHAO)
    omegadot    = -(GBEfficiency * Gen_Torque * gbRatio - ...
                        aerodynamicTorque)/(inertTot);
    WSE.omega   = WSE.omega + dt*omegadot;
    diff_omega  = - WSE.omega + Rotor_Speed * pi/30;
    WSE.Ee      = WSE.Ee + diff_omega * dt;
    WSE.V       = beta * WSE.Ee + gamma * diff_omega;
end
V_out = WSE.V;
end
