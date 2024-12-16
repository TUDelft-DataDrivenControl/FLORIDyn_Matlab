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

function VTKtoPNG(pathVTK,nameVTK,pathOutput,limits)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dirs = dir(pathVTK);

for iF = 3:length(dirs)
    path2VTK = [pathVTK '\' dirs(iF).name '\' nameVTK];
    
    % Load vtk data
    if strcmp(path2VTK(end-3:end),'.vtk')
        [~,cellCenters,cellData] = importVTK(path2VTK);
    else
        [~,cellCenters,cellData] = importVTK([path2VTK '.vtk']);
    end
    
    
    UmeanAbs = sqrt(sum(cellData.^2,2));
    X = cellCenters(:,1);
    Y = cellCenters(:,2);
    t=delaunay([X,Y]);
    
    % Plot result
    fig = figure(10);
    trisurf(t,X,Y,UmeanAbs,'EdgeColor','none')
    view(2)
    axis equal;
    axis tight;
    colormap(viridis(1000))
    caxis(limits)
    fig.Position = [29 319 1851 646];
    saveas(fig,[pathOutput '\' dirs(iF).name '.png'])
    
end
end

