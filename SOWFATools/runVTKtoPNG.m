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
pathVTK = 'W:\OpenFOAM\marcusbecker-2.4.0\simulationCases\2021_Paper_Marcus\3turb_yaw_lowTI\postProcessing\postProcessing\sliceDataInstantaneous';
pathOutput = 'C:\Users\marcusbecker\surfdrive\PhD_Surf\01_Research\01_FLORIDyn\02_Matlab\FLORIDynDDC\Simulations\2021_3T_Case_yaw_IandI\Data\HHFlowField';
nameVTK = 'U_slice_horizontal.vtk';
VTKtoPNG(pathVTK,nameVTK,pathOutput,[2,11])
pathOutput = 'C:\Users\marcusbecker\surfdrive\PhD_Surf\01_Research\01_FLORIDyn\02_Matlab\FLORIDynDDC\Simulations\2021_3T_Case_yaw_IandI\Data\HHFlowFieldAveraged';
VTKtoAvePNG(pathVTK,nameVTK,pathOutput,[3.5,9.5])