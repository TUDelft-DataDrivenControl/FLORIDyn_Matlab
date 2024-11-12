function VTKtoAvePNG(pathVTK,nameVTK,pathOutput,limits)
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
    
    % Horizontal slice trough the wake field (xy plane)
    % Create 2D interpolant
    interpolant = ...
        scatteredInterpolant(cellCenters(:,1:2),UmeanAbs);
    
    % x axis plot = x axis field
    Xaxis = linspace(...
        min(cellCenters(:,1),[],1),...
        max(cellCenters(:,1),[],1),900);
    
    % y axis plot = y axis field
    Yaxis = linspace(...
        min(cellCenters(:,2),[],1),...
        max(cellCenters(:,2),[],1),300);
    
    % Create meshgrid for interpolation
    [Xm,Ym] = meshgrid(Xaxis,Yaxis);
    UmeanAbs = interpolant(Xm,Ym);
    
    %I = cat(3,Xm,Ym,UmeanAbs);
    H_a = fspecial('average',[53,1]);
    
    If = imfilter(UmeanAbs,H_a);
    
    % Plot result
    fig = figure(12);
    %trisurf(t,X,Y,UmeanAbs,'EdgeColor','none')
    imagesc(Xaxis,Yaxis,If);
    set(gca,'YDir','normal');
    axis equal;
    axis tight;
    colormap(viridis(1000))
    caxis(limits)
    fig.Position = [29 319 1851 646];
    saveas(fig,[pathOutput '\' dirs(iF).name '.png'])
    
end
end

