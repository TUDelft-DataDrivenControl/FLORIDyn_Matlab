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

function paramFLORIS = parameterFLORIS()
%PARAMETERFLORIS returns all parameters needed to evaluate the FLORIS model
Bart            = false;
param_clc_2024  = false;
if Bart
    paramFLORIS.alpha   = 1.088; %#ok<UNRCH>
    paramFLORIS.beta    = 0.222;
    paramFLORIS.k_a     = 0.537;
    paramFLORIS.k_b     = -0.000848;
    paramFLORIS.k_fa    = 7.84;
    paramFLORIS.k_fb    = 4.57;
    paramFLORIS.k_fc    = 0.43;
    paramFLORIS.k_fd    = -0.246;
    % Note 2024/03/04: SOWFA power is off by a lot, so eta is higher than usual, otherwise too high wind speed. 
    paramFLORIS.eta     = 1.0%1.136; %1.39; %0.7918; %0.8572;
    paramFLORIS.p_p     = 2.2;
    paramFLORIS.airDen  = 1.225;    % Air density kg/m^3 (SOWFA)
    paramFLORIS.TIexp   = 3;        % TI crosswind distribution is modelled as 
                                    %   a factor of sigma_y and sigma_z. This
                                    %   is in contrast to previous literature
                                    %   which assumed a fixed box area and the
                                    %   turbines therein as source of higher TI
elseif param_clc_2024
    % Best sample
    %   FLORIDyn    0.6966    0.3570  206.2331    0.7396    
    %   FLORIS      pp 2.5433  alpha  0.7515  beta  0.3483  ka  0.7701  kfb  0.0467  kfc  0.2055  kTI  1.5374    
    %   EnKF        6.8011    0.1991    0.2633
    paramFLORIS.alpha   = 0.7515;
    paramFLORIS.beta    = 0.3483;
    paramFLORIS.k_a     = 0.7701;
    paramFLORIS.k_b     = 0.003678;
    paramFLORIS.k_fa    = 0.73;
    paramFLORIS.k_fb    = 0.0467;
    paramFLORIS.k_fc    = 0.2055;
    paramFLORIS.k_fd    = -0.32;
    % Note 2024/03/04: SOWFA power is off by a lot, so eta is higher than usual, otherwise too high wind speed. 
    paramFLORIS.eta     = 1.136; %1.39; %0.7918; %0.8572;
    paramFLORIS.p_p     = 2.5433;
    paramFLORIS.airDen  = 1.225;    % Air density kg/m^3 (SOWFA)
    paramFLORIS.TIexp   = 1.5374;        % TI crosswind distribution is modelled as 
                                    %   a factor of sigma_y and sigma_z. This
                                    %   is in contrast to previous literature
                                    %   which assumed a fixed box area and the
                                    %   turbines therein as source of higher TI
else
    paramFLORIS.alpha   = 2.32;
    paramFLORIS.beta    = 0.154;
    paramFLORIS.k_a     = 0.38371;
    paramFLORIS.k_b     = 0.003678;
    paramFLORIS.k_fa    = 0.73;
    paramFLORIS.k_fb    = 0.8325;
    paramFLORIS.k_fc    = 0.0325;
    paramFLORIS.k_fd    = -0.32;
    % Note 2024/03/04: SOWFA power is off by a lot, so eta is higher than usual, otherwise too high wind speed. 
    paramFLORIS.eta     = 1.0; %1.39; %0.7918; %0.8572;
    paramFLORIS.p_p     = 2.2;
    paramFLORIS.airDen  = 1.225;    % Air density kg/m^3 (SOWFA)
    paramFLORIS.TIexp   = 3;        % TI crosswind distribution is modelled as 
                                    %   a factor of sigma_y and sigma_z. This
                                    %   is in contrast to previous literature
                                    %   which assumed a fixed box area and the
                                    %   turbines therein as source of higher TI
end


end

