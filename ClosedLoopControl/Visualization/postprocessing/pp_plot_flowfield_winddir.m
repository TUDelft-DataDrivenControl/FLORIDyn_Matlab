function pp_plot_flowfield_winddir(path)
%PP_PLOT_FLOWFIELD_WINDDIR Summary of this function goes here
%   Detailed explanation goes here

% Hardcoded settings from setup.mlx
nOP                  = 200;
Sigma.IterSigma_DW   = 512;%512;
Sigma.IterSigma_CW   = 512; %63
Sigma.IterSigma_time = 50;%126;
timeStep             = 5;

data_dirs = dir([path filesep 'FlowField']);

mkdir([path filesep 'FlowField_plots2'])

X = readmatrix([path filesep 'FlowField' filesep 'X.csv']);
Y = readmatrix([path filesep 'FlowField' filesep 'Y.csv']);

[X_Wind, Y_Wind] = meshgrid( ...
    linspace(X(1,1),X(1,end),20), ...
    linspace(Y(1,1),Y(end,1),20));

% X_Wind(:,1) = []; X_Wind(:,end) = [];
% Y_Wind(:,1) = []; Y_Wind(:,end) = [];
% X_Wind(1,:) = []; X_Wind(end,:) = [];
% Y_Wind(1,:) = []; Y_Wind(end,:) = [];

OPx = readmatrix([path filesep 'FlowField' filesep 'OP_x.csv']);
OPy = readmatrix([path filesep 'FlowField' filesep 'OP_y.csv']);

dir_states = readmatrix([path filesep 'EnKF_States' filesep 'dir.txt']);
dir_data   = readmatrix([path filesep '..' filesep '..' filesep 'LESData/dir.csv']);
dir_data(:,1) = dir_data(:,1) + 30000;
%% Get bounds on the colormap
u_max = 0;
u_min = 1000;

for iD = 3:2:length(data_dirs)-4
    data_U = readmatrix([data_dirs(iD).folder filesep data_dirs(iD).name]);
    u_max = max(u_max, max(data_U,[],'all'));
    u_min = min(u_min, min(data_U,[],'all'));
end



%
i = 0;
for iD = 3:2:length(data_dirs)-4
    i = i+1;
    data_U = readmatrix([data_dirs(iD).folder filesep data_dirs(iD).name]);
    data_W = readmatrix([data_dirs(iD+1).folder filesep data_dirs(iD+1).name]);


    W = getWeightsPar( ...
        X_Wind(:), Y_Wind(:), OPx(i,2:end)', OPy(i,2:end)', ...
        nOP, 10, Sigma, timeStep);
    W = W./sum(W,2);

    phi_grid = reshape(W*dir_states(i,:)',size(X_Wind));
    phi_grid = angSOWFA2world(phi_grid);
    phi_true = angSOWFA2world(interp1(dir_data(:,1), dir_data(:,2), OPx(i,1)));

    figure
    %contourf(X,Y,data_U,linspace(u_min,u_max,30),'EdgeColor','none')
    imagesc(X(1,:)',Y(:,1),data_U)
    colormap((Bu(30)))
    hold on
    % EnKF_plot_mask(X,Y,-data_W,-.2, [1,1,1])
    scatter(OPx(i,2:end),OPy(i,2:end),2,"white",'filled')
    scatter(OPx(i,2:10:end),OPy(i,2:10:end),10,"white",'filled')
    quiver(X_Wind,Y_Wind, ...
        ones(size(X_Wind))*cos(phi_true), ...
        ones(size(X_Wind))*sin(phi_true), ...
        .71,'filled','Color','w')
    quiver(X_Wind,Y_Wind,cos(phi_grid), sin(phi_grid),.71,'filled','Color',[0.0196    0.1882    0.3804])
    hold off
    axis equal
    xlabel('Easting in m')
    ylabel('Northing in m')
    c = colorbar;
    c.Label.String = 'Wind speed in m/s';
    title(['Time: ' data_dirs(iD).name(1:end-6)])
    clim([u_min,u_max])
    xlim([X(1,1),X(1,end)])
    ylim([Y(1,1),Y(end,1)])
    set(gca,'YDir','normal')
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'Masked_FlowField' data_dirs(iD).name(1:end-6) '.jpg'],...
        'Resolution',200)
    % exportgraphics(gcf,...
    %     [path filesep 'FlowField_plots' filesep ...
    %     'Masked_FlowField' data_dirs(iD).name(1:end-6) '.pdf'],...
    %     'ContentType','vector')
    close(gcf)


    % figure
    % contourf(X,Y,data_U,linspace(u_min,u_max,30),'EdgeColor','none')
    % colormap((Bu(30)))
    % hold on
    % scatter(OPx(i,2:end),OPy(i,2:end),5,"white",'filled')
    % EnKF_plot_mask(X,Y,-data_W,-.2, [1,1,1])
    % %quiver(X_Wind,Y_Wind,cos(phi_grid), sin(phi_grid),.71,'filled','Color','w')
    % hold off
    % axis equal
    % xlabel('Easting in m')
    % ylabel('Northing in m')
    % c = colorbar;
    % c.Label.String = 'Wind speed in m/s';
    % xlim([X(1,1), X(1,end)])
    % ylim([Y(1,1), Y(end,1)])
    % title(['Time: ' data_dirs(iD).name(1:end-6)])
    % clim([u_min,u_max])
    % xlim([X(1,1),X(1,end)])
    % ylim([Y(1,1),Y(end,1)])
    % 
    % 
    % exportgraphics(gcf,...
    %     [path filesep 'FlowField_plots' filesep ...
    %     'Masked_FlowField_w_OPs' data_dirs(iD).name(1:end-6) '.jpg'],...
    %     'Resolution',200)
    % exportgraphics(gcf,...
    %     [path filesep 'FlowField_plots' filesep ...
    %     'Masked_FlowField_w_OPs' data_dirs(iD).name(1:end-6) '.pdf'],...
    %     'ContentType','vector')
    % close(gcf)
end





end

function W = getWeightsPar(GP_x, GP_y, States_OP_x, States_OP_y, nOP, nT, Sigma, timeStep)
% GETWEIGHTS generates a matrix with weights based on the distance of the
% OPs to each other. Wind direction is taken into account, also a temporal
% decay factor is implemented. Implementation inspired by Lejeune.

dist_OP_2_GP = sqrt((GP_x - States_OP_x').^2 + (GP_y - States_OP_y').^2);
W = exp(dist_OP_2_GP./(2 * Sigma.IterSigma_DW^2));

W = W .* repmat(...
    exp(-((0:nOP-1)*timeStep).^2 / (2 * Sigma.IterSigma_time^2)),...
    1,nT);
end