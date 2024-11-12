% Postprocess Data
path = '/Users/marcusbecker/surfdrive/PhD_Surf/01_Research/01_FLORIDyn/02_Matlab/marcusbecker@hpc06.tudelft.net/Simulations/2024_CLC/10T_HKN_het_p10ramp_201deg/Results/first_run';

%% Plot trajectories
close all
wind_dir = readmatrix([path filesep 'EnKF_States' filesep 'dir.txt']);

% REPLACE!!
enkf_t = 30030:30:31980;

data_dirs = dir([path filesep 'Con_Trajectories']);
nT        = 0;
t_min     = 0;
t_max     = 0;
delta_t   = 0;
yaw_min   = 1000;
yaw_max   = -1000;
ah        = str2double(data_dirs(4).name(1:end-4)) - ...
    str2double(data_dirs(3).name(1:end-4));

cols    = Bu(4);
col     = cols(3,:);
cols    = RdBu(7);

for iD = 3:length(data_dirs)
    data = readmatrix([data_dirs(iD).folder filesep data_dirs(iD).name]);
    if nT == 0; nT = size(data,2)-1; end
    if delta_t == 0
        delta_t = data(2,1)-data(1,1); 
        ah      = ah/delta_t;
    end
    if iD == 3; t_min = data(1,1); end
    if iD == length(data_dirs); t_max = data(end,1); end
    
    yaw_max = max(yaw_max,max(data(:,2:end),[],'all'));
    yaw_min = min(yaw_min,min(data(:,2:end),[],'all'));

    for iT = 1:nT
        figure(iT)
        hold on
        plot(data(1:end,1)/delta_t,data(:,iT+1),...
            "Color",[.7,.7,.7],'LineStyle','-','LineWidth',2)
        plot(data(1:ah,1)/delta_t,data(1:ah,iT+1),...
            "Color",'k','LineStyle','-','LineWidth',2)
    end
end

mkdir([path filesep 'Con_Trajectories_plots'])
nOP = size(wind_dir,2)/nT;
for iT = 1:nT
    figure(iT)
    hold on
    yaw_ref = 30:-10:-30;
    for iy = 1:length(yaw_ref)
        plot(enkf_t/delta_t,yaw_ref(iy)+wind_dir(:,(iT-1)*nOP+1),...
            'Color',cols(iy,:),'LineStyle','--','LineWidth',1.5)
    end
    hold off
    title(['Turbine ' num2str(iT)])
    xlabel('Time Step')
    ylabel('Orientation')
    xlim([t_min, t_max]/delta_t)
    ylim([yaw_min, yaw_max])
    grid on

    % Export
    % exportgraphics(gcf,...
    %     [path filesep 'Con_Trajectories_plots' filesep ...
    %     'T' pad(num2str(iT),2,'left','0') '.jpg'],...
    %     'Resolution',200)
    % exportgraphics(gcf,...
    %     [path filesep 'Con_Trajectories_plots' filesep ...
    %     'T' pad(num2str(iT),2,'left','0') '.pdf'],...
    %     'ContentType','vector')
    close(gcf)
end

%% Plot landscape
data_dirs = dir([path filesep 'Con_Landscape']);
mkdir([path filesep 'Con_Landscape_plots'])

for iD = 3:3:length(data_dirs)
    E1 = readmatrix([data_dirs(iD).folder filesep data_dirs(iD).name]);
    E2 = readmatrix([data_dirs(iD+1).folder filesep data_dirs(iD+1).name]);
    J  = readmatrix([data_dirs(iD+2).folder filesep data_dirs(iD+2).name]);

    J(J==0) = nan;

    figure
    contourf(E1,E2,J,100,'EdgeColor','none')
    hold on
    scatter(0,0,40,'white')
    axis equal
    xlabel('e_1')
    ylabel('e_2')
    title(['Cost function, t=' data_dirs(iD).name(1:end-7) ' s'])
    colormap(flipud(Bu(100)))
    grid on

    % Export
    exportgraphics(gcf,...
        [path filesep 'Con_Landscape_plots' filesep ...
        'Landscape' data_dirs(iD).name(1:end-7) '.jpg'],...
        'Resolution',200)
    exportgraphics(gcf,...
        [path filesep 'Con_Landscape_plots' filesep ...
        'Landscape' data_dirs(iD).name(1:end-7) '.pdf'],...
        'ContentType','vector')
    close(gcf)
end

%% Get and plot time statistics
fix = true; % ONLY NEEDED FOR THE FIRST RUN (7th of May 2024)

data_dirs = dir([path filesep 'Con_Optimization_groups']);
mkdir([path filesep 'Time_plots'])

time_per_turbine = [];
time_total = [];
for iD = 3:length(data_dirs)
    data = readmatrix([data_dirs(iD).folder filesep data_dirs(iD).name]);
    
    if fix
        data(2:end,end) = data(2:end,end) - ...
            data(1:end-1,end);
    end
    
    time_per_turbine = [time_per_turbine;...
        sum(data(:,1:end-1),2), data(:,end)./sum(data(:,1:end-1),2)];

    time_total = [time_total;...
        sum(data(:,1:end-1),2), data(:,end)];
end
%
figure
boxchart(time_total(:,1),time_total(:,2))
xlabel('Number of turbines in Optimization')
ylabel('Time in s')
xticks(nTpO)
grid on

% Export
exportgraphics(gcf,...
    [path filesep 'Time_plots' filesep ...
    'Optimization_time.jpg'],...
    'Resolution',200)
exportgraphics(gcf,...
    [path filesep 'Time_plots' filesep ...
    'Optimization_time.pdf'],...
    'ContentType','vector')
close(gcf)

%%
nTpO = unique(time_per_turbine(:,1));
figure
boxchart(time_per_turbine(:,1),time_per_turbine(:,2))
xlabel('Number of turbines in Optimization')
ylabel('Time per turbine in s')
xticks(nTpO)
grid on
% Export
exportgraphics(gcf,...
    [path filesep 'Time_plots' filesep ...
    'Optimization_time_per_turbine.jpg'],...
    'Resolution',200)
exportgraphics(gcf,...
    [path filesep 'Time_plots' filesep ...
    'Optimization_time_per_turbine.pdf'],...
    'ContentType','vector')
close(gcf)

%% Norm by 1T case (if present)
if nTpO(1) == 1
    mean_1T = mean(time_per_turbine(time_per_turbine(:,1)==1,2));
    figure
    for i = 1:length(nTpO)
        boxchart(time_per_turbine(:,1),time_per_turbine(:,2)/mean_1T)
    end
    grid on
    xlabel('Number of turbines in Optimization')
    ylabel('Time per turbine in s / mean time for one turbine')
    xticks(nTpO)

    exportgraphics(gcf,...
        [path filesep 'Time_plots' filesep ...
        'Optimization_time_per_turbine_norm.jpg'],...
        'Resolution',200)
    exportgraphics(gcf,...
        [path filesep 'Time_plots' filesep ...
        'Optimization_time_per_turbine_norm.pdf'],...
        'ContentType','vector')
    close(gcf)
end

%%
data = readmatrix([path filesep 'EnKF_Time' filesep 'time_sim_com_cor.csv']);

figure
boxchart(data(2:end,2:end))
xlabel('EnKF steps')
ylabel('Time in s')
xticklabels({'Simulation','Combination','Correction'})

exportgraphics(gcf,...
    [path filesep 'Time_plots' filesep ...
    'EnKF_Time.jpg'],...
    'Resolution',200)
exportgraphics(gcf,...
    [path filesep 'Time_plots' filesep ...
    'EnKF_Time.pdf'],...
    'ContentType','vector')
close(gcf)

%% Plot flow fields
data_dirs = dir([path filesep 'FlowField']);

mkdir([path filesep 'FlowField_plots'])

X = readmatrix([path filesep 'FlowField' filesep 'X.csv']);
Y = readmatrix([path filesep 'FlowField' filesep 'Y.csv']);

OPx = readmatrix([path filesep 'FlowField' filesep 'OP_x.csv']);
OPy = readmatrix([path filesep 'FlowField' filesep 'OP_y.csv']);
i = 0;
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

    figure
    contourf(X,Y,data_U,linspace(u_min,u_max,30),'EdgeColor','none')
    colormap((Bu(30)))
    hold on
    EnKF_plot_mask(X,Y,-data_W,-.2, [1,1,1])
    hold off
    axis equal
    xlabel('Easting in m')
    ylabel('Northing in m')
    c = colorbar;
    c.Label.String = 'Wind speed in m/s';
    title(['Time: ' data_dirs(iD).name(1:end-6)])

    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'Masked_FlowField.jpg'],...
        'Resolution',200)
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'Masked_FlowField.pdf'],...
        'ContentType','vector')
    close(gcf)


    figure
    contourf(X,Y,data_U,linspace(u_min,u_max,30),'EdgeColor','none')
    colormap((Bu(30)))
    hold on
    scatter(OPx(i,2:end),OPy(i,2:end),5,"white",'filled')
    EnKF_plot_mask(X,Y,-data_W,-.2, [1,1,1])
    hold off
    axis equal
    xlabel('Easting in m')
    ylabel('Northing in m')
    c = colorbar;
    c.Label.String = 'Wind speed in m/s';
    xlim([X(1,1), X(1,end)])
    ylim([Y(1,1), Y(end,1)])
    title(['Time: ' data_dirs(iD).name(1:end-6)])
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'Masked_FlowField_w_OPs.jpg'],...
        'Resolution',200)
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'Masked_FlowField_w_OPs.pdf'],...
        'ContentType','vector')
    close(gcf)
end
%% Plot OP traces

figure
hold on
for i = 2:size(OPx,2)
    plot(OPx(:,i),OPy(:,i),'-k')
end
axis equal
xlim([X(1,1), X(1,end)])
ylim([Y(1,1), Y(end,1)])
xticks([])
yticks([])
title('OP positions')
exportgraphics(gcf,...
    [path filesep 'FlowField_plots' filesep ...
    'OP_Positions.jpg'],...
    'Resolution',200)
exportgraphics(gcf,...
    [path filesep 'FlowField_plots' filesep ...
    'OP_Positions.pdf'],...
    'ContentType','vector')
close(gcf)

nOP = (size(OPx,2)-1)/nT;
for iT = 1:nT
    figure
    hold on
    for iOP = 1:nOP
        plot(OPx(:,iOP+(iT-1)*nOP),OPy(:,iOP+(iT-1)*nOP),'-k')
    end
    axis equal
    xlim([X(1,1), X(1,end)])
    ylim([Y(1,1), Y(end,1)])
    xticks([])
    yticks([])
    title(['OP positions Turbine ' num2str(iT)])
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'OP_Positions_T' num2str(iT) '.jpg'],...
        'Resolution',200)
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'OP_Positions_T' num2str(iT) '.pdf'],...
        'ContentType','vector')
    close(gcf)
end
