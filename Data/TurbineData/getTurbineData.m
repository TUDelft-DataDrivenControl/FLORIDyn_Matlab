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


function t_data = getTurbineData(name)
%GETTURBINEDATA stores the turbine type specific information. It expects a
%one dimensional cell struct with strings and will return a struct with the
%according data (Hubheight, Diameter)
t_data.h  = zeros(length(name),3);
t_data.D  = zeros(length(name),1);

for i = 1:length(name)
    switch name{i}
        case 'DTU 10MW'
            % Nacelle position seen from base
            t_data.NacPos(i,:) = [0,0,119];
            % Diameter
            t_data.D(i) = 178.4;
        case 'DTU 5MW'
            % Hubheight
            t_data.NacPos(i,:) = [0,0,119];      % REPLACE
            % Diameter
            t_data.D(i) = 178.4;    % REPLACE
        case 'Senvion 6.2M'
            % Sources
            %   https://en.wind-turbine-models.com/turbines/885-senvion-6.2m126-offshore
            %   https://www.nordseeone.com/engineering-construction/wind-turbines/general.html
            % Nacelle position seen from base
            t_data.NacPos(i,:) = [0,0,152-29];
            % Diameter
            t_data.D(i) = 126;
        case 'V116'
            % Source
            %   https://en.wind-turbine-models.com/turbines/1814-vestas-v116-2.1
            % Nacelle position seen from base
            t_data.NacPos(i,:) = [0,0,84]; %or 91.5
            % Diameter
            t_data.D(i) = 116;
        case 'V117'
            % Source
            %   https://www.vestas.com/en/products/4-mw-platform/V117-4-2-MW
            % Nacelle position seen from base
            t_data.NacPos(i,:) = [0,0,84]; %or 91.5
            % Diameter
            t_data.D(i) = 117;
        case 'V162'
            % Source
            %   https://www.vestas.com/en/products/enventus-platform/v162-6-8-mw
            % Nacelle position seen from base
            t_data.NacPos(i,:) = [0,0,119]; %or 169
            % Diameter
            t_data.D(i) = 162;
        case 'GE Haliade X'
            % Source
            %   https://www.thewindpower.net/turbine_en_1579_ge-energy_haliade-x-12-mw.php
            % Nacelle position seen from base
            t_data.NacPos(i,:) = [0,0,150];
            % Diameter
            t_data.D(i) = 220;
            
        otherwise
            error(['Turbine type ' name{i} ' not known or misspelled.'])
    end
end
end

