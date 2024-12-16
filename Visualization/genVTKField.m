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

function genVTKField(x,y,z,m,fileName,metricName)
%genVTKField creates a vtk file from 3d data which can be visualized with
%Paraview for instance.

if strcmp(fileName(end-4:end),'.vtk')
    filename=fileName;
else
    filename=[fileName '.vtk'];
end

nr_of_elements=numel(x);
fid = fopen(filename, 'w'); 
%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'BINARY\n\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(x,2)) ' ' num2str(size(x,3)) '\n']);
fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
fclose(fid);
%append binary x,y,z data
fid = fopen(filename, 'a'); 
fwrite(fid, [reshape(x,1,[]);  reshape(y,1,[]); reshape(z,1,[])],'float','b');
%append another ASCII sub header
fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);

%append measurement data
fprintf(fid, ['\nSCALARS ' metricName ' float\n']); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fwrite (fid, reshape(m,1,[]),'float','b'); %binary data
fclose(fid);

end

