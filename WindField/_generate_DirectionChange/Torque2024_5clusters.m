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

function Torque2024_5clusters(dir,t_0,t_start,t_end,path)
%TORQUE2024_5CLUSTERS Generates the five clusters wind direction identified
%for the Torque2024 conference

cluster_o_deg = [
    0, 0.28, 0.71, 1.16, 1.62, 2.07, 2.51, 2.94, 3.34, 3.7,  4.02;
    0, 0.11, 0.27, 0.44, 0.62, 0.8,  0.97, 1.13, 1.29, 1.44, 1.57;
    0, 0,    0,    0,    0,    0,    0,    0,    0,    0,    0;
    0, -.11, -.29, -.47, -.66, -.84, -1.02, -1.19, -1.34, -1.49, -1.61;
    0, -0.28, -0.71, -1.18, -1.64, -2.1, -2.55, -2.99, -3.39, -3.77, -4.1];
cluster_t = [0, (1:10) * 60];

%%
if ~endsWith(path,filesep)
    path = [path filesep];
end
%% Generate clusters
for iC = 1:size(cluster_o_deg,1)
    winddir = [cluster_t'+t_start cluster_o_deg(iC,:)'+dir];
    winddir = [winddir(1,:);winddir;winddir(end,:)];
    winddir(1,1)    = t_0;
    winddir(1,end)  = t_end;

    writematrix(winddir,[path 'WindDirC' num2str(iC) '_' num2str(dir) ...
        '.csv'])
end


end

