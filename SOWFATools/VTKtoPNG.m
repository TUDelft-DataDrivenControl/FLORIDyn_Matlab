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

