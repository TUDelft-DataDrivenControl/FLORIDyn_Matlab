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

pathToSimulation = ['2024_CLC' filesep '10T_HKN_LUT'];
Dirs             = 1:2:360;

lut = createLUTwFLORIDyn(Dirs, pathToSimulation);

writematrix(lut,'LuT_orientation_HKN_no_rotation.csv')

%%
lut_ori = readmatrix("Simulations/2024_CLC/10T_HKN_LUT/Results/HKN_LuT_orientation.csv");

lut_yaw = Dirs' - lut_ori;
deg_yaw_table = [Dirs' lut_yaw];
deg_ori_table = [Dirs' lut_ori];

%%
writematrix(deg_ori_table,'/Users/marcusbecker/surfdrive/PhD_Surf/02_Communication/03_Papers/07_CL_Control/Reference Controllers/Yaw steering/deg_ori_table_internal.csv')
writematrix(deg_yaw_table,'/Users/marcusbecker/surfdrive/PhD_Surf/02_Communication/03_Papers/07_CL_Control/Reference Controllers/Yaw steering/deg_yaw_table_internal.csv')
%%
figure
imagesc(lut_yaw)
colormap(flip(RdBu(61)))
clim([-30.5,30.5])