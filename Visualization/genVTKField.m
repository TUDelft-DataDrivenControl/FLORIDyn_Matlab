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

