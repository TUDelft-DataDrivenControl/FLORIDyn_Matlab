function postprocess_output(path)
%POSTPROCESS_OUTPUT Summary of this function goes here
%   Detailed explanation goes here

% path = '/marcusbecker@hpc06.tudelft.net/Simulations/2024_CLC/10T_HKN_het_p10ramp_201deg/Results/first_run';
% try
    pp_plot_flowfield_winddir(path)
    nT = private_trajectories(path);
    % private_landscape(path);
    private_time_stats(path);
    private_flow_fields(path);
    private_plot_OP_traces(path,nT);
    private_enkf_states(path, nT)
% catch
%     disp('Error during postprocessing!!')
% end

end



function private_plot_OP_traces(path,nT)

X = readmatrix([path filesep 'FlowField' filesep 'X.csv']);
Y = readmatrix([path filesep 'FlowField' filesep 'Y.csv']);

OPx = readmatrix([path filesep 'FlowField' filesep 'OP_x.csv']);
OPy = readmatrix([path filesep 'FlowField' filesep 'OP_y.csv']);

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
    title(['OP positions Turbine ' num2str(iT-1)])
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'OP_Positions_T' num2str(iT-1) '.jpg'],...
        'Resolution',200)
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'OP_Positions_T' num2str(iT-1) '.pdf'],...
        'ContentType','vector')
    close(gcf)
end
end

function private_flow_fields(path)

data_dirs = dir([path filesep 'FlowField']);

mkdir([path filesep 'FlowField_plots'])

X = readmatrix([path filesep 'FlowField' filesep 'X.csv']);
Y = readmatrix([path filesep 'FlowField' filesep 'Y.csv']);

OPx = readmatrix([path filesep 'FlowField' filesep 'OP_x.csv']);
OPy = readmatrix([path filesep 'FlowField' filesep 'OP_y.csv']);
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
    clim([u_min,u_max])

    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'Masked_FlowField' data_dirs(iD).name(1:end-6) '.jpg'],...
        'Resolution',200)
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'Masked_FlowField' data_dirs(iD).name(1:end-6) '.pdf'],...
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
    clim([u_min,u_max])
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'Masked_FlowField_w_OPs' data_dirs(iD).name(1:end-6) '.jpg'],...
        'Resolution',200)
    exportgraphics(gcf,...
        [path filesep 'FlowField_plots' filesep ...
        'Masked_FlowField_w_OPs' data_dirs(iD).name(1:end-6) '.pdf'],...
        'ContentType','vector')
    close(gcf)
end

end

function private_time_stats(path)
cols = Bu(4);

data_dirs = dir([path filesep 'Con_Optimization_groups']);
mkdir([path filesep 'Time_plots'])

time_per_turbine = [];
time_total = [];
for iD = 3:length(data_dirs)
    data = readmatrix([data_dirs(iD).folder filesep data_dirs(iD).name]);
    
    time_per_turbine = [time_per_turbine;...
        sum(data(:,1:end-1),2), data(:,end)./sum(data(:,1:end-1),2)];

    time_total = [time_total;...
        sum(data(:,1:end-1),2), data(:,end)];
end
%
figure
boxchart(time_total(:,1),time_total(:,2), ...
    "BoxFaceColor",cols(2,:),'MarkerColor',cols(2,:),'BoxEdgeColor',cols(2,:))
xlabel('Number of turbines in Optimization')
ylabel('Time in s')
nTpO = unique(time_per_turbine(:,1));
xticks(nTpO)
grid on
title('Total optimization time in s')

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

%

figure
boxchart(time_per_turbine(:,1),time_per_turbine(:,2), ...
    "BoxFaceColor",cols(3,:),'MarkerColor',cols(3,:),'BoxEdgeColor',cols(3,:))
xlabel('Number of turbines in Optimization')
ylabel('Time per turbine in s')
xticks(nTpO)
grid on
title('Optimization time per turbine in s')
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

% Norm by 1T case (if present)
if nTpO(1) == 1
    mean_1T = mean(time_per_turbine(time_per_turbine(:,1)==1,2));
    figure
    for i = 1:length(nTpO)
        boxchart(time_per_turbine(:,1),time_per_turbine(:,2)/mean_1T, ...
            "BoxFaceColor",cols(4,:),'MarkerColor',cols(4,:),'BoxEdgeColor',cols(4,:))
    end
    grid on
    xlabel('Number of turbines in Optimization')
    ylabel('Time per turbine in s / mean time for one turbine')
    xticks(nTpO)
    title('Optimization time normalized by one turbine time in s')
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

% ===== EnKF time
data = readmatrix([path filesep 'EnKF_Time' filesep 'time_sim_com_cor.csv']);
cols = Rd(2);
figure
boxchart(data(:,2:end),...
    'BoxFaceColor',cols(1,:),'MarkerColor',cols(1,:),'BoxEdgeColor',cols(1,:))
title('EnKF time statistics')
xticklabels({'Simulation', 'Combination', 'Correction'})
ylabel('Time in s')
grid on

exportgraphics(gcf,...
    [path filesep 'Time_plots' filesep ...
    'EnKF_time.jpg'],...
    'Resolution',200)
exportgraphics(gcf,...
    [path filesep 'Time_plots' filesep ...
    'EnKF_time.pdf'],...
    'ContentType','vector')
close(gcf)

end

function private_landscape(path)

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

end

function nT = private_trajectories(path)

close all
wind_dir    = readmatrix([path filesep 'EnKF_States' filesep 'dir.txt']);
enkf_time   = readmatrix([path filesep 'EnKF_Time' filesep 'time_sim_com_cor.csv']);
enkf_t      = enkf_time(:,1);
enkf_t      = enkf_t(and(enkf_t>31200, mod(enkf_t,60) == 0));

clear enkf_time
yaw_SOWFA = readmatrix([path filesep 'turbineOutput' filesep 'nacYaw']);

data_dirs = dir([path filesep 'Con_Trajectories']);
nT        = 0;
t_min     = 0;
t_max     = 0;
delta_t   = 0;
yaw_min   = min(yaw_SOWFA(:,4));
yaw_max   = max(yaw_SOWFA(:,4));
ah        = str2double(data_dirs(4).name(1:end-4)) - ...
    str2double(data_dirs(3).name(1:end-4));
nOP       = 0;
cols    = RdBu(7);
cols(4,:) = [.5 .5 .5];
yaw_ref = 30:-10:-30;

for iD = 3:length(data_dirs)
    data = readmatrix([data_dirs(iD).folder filesep data_dirs(iD).name]);
    if nT == 0; nT = size(data,2)-1; end
    if delta_t == 0
        delta_t = data(2,1)-data(1,1); 
        ah      = ah/delta_t;
    end
    if iD == 3; t_min = data(1,1); end
    if iD == length(data_dirs); t_max = data(end,1); end
    if nOP == 0; nOP = size(wind_dir,2)/nT; end


    yaw_max = max(yaw_max,max(data(:,2:end),[],'all'));
    yaw_min = min(yaw_min,min(data(:,2:end),[],'all'));

    for iT = 1:nT
        figure(iT)
        hold on
        for iy = 1:length(yaw_ref)
            plot(enkf_t(end-size(wind_dir,1)+1:end),yaw_ref(iy)+wind_dir(:,(iT-1)*nOP+1),...
                'Color',cols(iy,:),'LineStyle','-','LineWidth',1.5)
        end
        plot(data(1:end,1),data(:,iT+1),...
            "Color",[.8,.8,.8],'LineStyle','-','LineWidth',1)
        plot(data(1:ah+1,1),data(1:ah+1,iT+1),...
            "Color",'k','LineStyle','-','LineWidth',2)
    end
end

mkdir([path filesep 'Con_Trajectories_plots'])

for iT = 1:nT
    figure(iT)
    hold on
    plot(yaw_SOWFA(iT:nT:end,2), yaw_SOWFA(iT:nT:end,4),...
        'LineWidth',1,'Color',cols(2,:))
    hold off
    title(['Turbine ' num2str(iT-1)])
    xlabel('Time Step')
    ylabel('Orientation')
    xlim([t_min, t_max])
    ylim([yaw_min-5, yaw_max+5])
    grid on

    % Export
    exportgraphics(gcf,...
        [path filesep 'Con_Trajectories_plots' filesep ...
        'T' pad(num2str(iT-1),2,'left','0') '.jpg'],...
        'Resolution',200)
    exportgraphics(gcf,...
        [path filesep 'Con_Trajectories_plots' filesep ...
        'T' pad(num2str(iT-1),2,'left','0') '.pdf'],...
        'ContentType','vector')
    close(gcf)
end

end


function private_enkf_states(path, nT)
%% Plot EnKF states
dir_data = ...
    readmatrix([path filesep 'EnKF_States' filesep 'dir.txt']);
dir_var_data = ...
    readmatrix([path filesep 'EnKF_States' filesep 'dir_var.txt']);
vel_data = ...
    readmatrix([path filesep 'EnKF_States' filesep 'vel.txt']);
vel_var_data = ...
    readmatrix([path filesep 'EnKF_States' filesep 'vel_var.txt']);

nOP = size(dir_data,2)/nT;
%% == Vel
figure
subplot(2,1,1)
imagesc(vel_data)
hold on
for iT = 1:nT-1
    plot([1,1]*nOP*iT+.5, [.5,size(vel_data,1)+.5],'-k','LineWidth',1)
end
xlabel('States per turbine')
ylabel('EnKF step')
xticks(nOP/2:nOP:nT*nOP)
xticklabels(0:nT-1)
colormap(Bu)
c = colorbar;
c.Label.String = 'Wind speed in m/s';
title('Estimated velocity')

subplot(2,1,2)
imagesc(vel_var_data)
hold on
for iT = 1:nT-1
    plot([1,1]*nOP*iT+.5, [.5,size(vel_data,1)+.5],'-k','LineWidth',1)
end
xlabel('States per turbine')
ylabel('EnKF step')
xticks(nOP/2:nOP:nT*nOP)
xticklabels(0:nT-1)
colormap(Bu)
c = colorbar;
c.Label.String = 'Variance in wind speed in m/s';
title('Variance of the estimated velocity')

mkdir([path filesep 'EnKF_States_Plots'])

exportgraphics(gcf,...
    [path filesep 'EnKF_States_Plots' filesep ...
    'Wind_speed.jpg'],...
    'Resolution',200)
exportgraphics(gcf,...
    [path filesep 'EnKF_States_Plots' filesep ...
    'Wind_speed.pdf'],...
    'ContentType','vector')
close(gcf)

%% == Dir
figure
subplot(2,1,1)
imagesc(dir_data)
hold on
for iT = 1:nT-1
    plot([1,1]*nOP*iT+.5, [.5,size(vel_data,1)+.5],'-k','LineWidth',1)
end
xlabel('States per turbine')
ylabel('EnKF step')
xticks(nOP/2:nOP:nT*nOP)
xticklabels(0:nT-1)
colormap(flipud(Rd))
c = colorbar;
c.Label.String = 'Wind direction in deg';
title('Estimated wind direction in deg')

subplot(2,1,2)
imagesc(dir_var_data)
hold on
for iT = 1:nT-1
    plot([1,1]*nOP*iT+.5, [.5,size(vel_data,1)+.5],'-k','LineWidth',1)
end
xlabel('States per turbine')
ylabel('EnKF step')
xticks(nOP/2:nOP:nT*nOP)
xticklabels(0:nT-1)
colormap(flipud(Rd))
c = colorbar;
c.Label.String = 'Variance in wind direction in deg';
title('Variance of the estimated wind direction')


exportgraphics(gcf,...
    [path filesep 'EnKF_States_Plots' filesep ...
    'Wind_direction.jpg'],...
    'Resolution',200)
exportgraphics(gcf,...
    [path filesep 'EnKF_States_Plots' filesep ...
    'Wind_direction.pdf'],...
    'ContentType','vector')

end

function private_EnKF_power(path, nT)
M = readtable([path filesep 'EnKF_Measurements' filesep 'EnKF_measurements.csv']);
powSOWFA = readmatrix([path filesep 'turbineOutput' filesep 'powerGenerator']);
powSOWFA(:,4) = powSOWFA(:,4)./(1.225*10^6);
%
maxP = max([M.PowerGenerated_MW_;powSOWFA(powSOWFA(:,4)>0,4)]);
minP = min([M.PowerGenerated_MW_;powSOWFA(powSOWFA(:,4)>0,4)]);

cols    = Bu(3);
nE      = sum(M.Time_s_==M.Time_s_(1))/nT;
nt      = length(unique(M.Time_s_));
p_mean = mean(reshape(M.PowerGenerated_MW_,[],nE),2);

for iT = 1:nT
    figure
    hold on
    for iE = 1:nE
        plot(M.Time_s_(iT:nT:nt*nT),...
            M.PowerGenerated_MW_((iE-1)*nt*nT+iT:nT:iE*nt*nT), ...
            'Color',cols(2,:),'LineWidth',.5)
    end
    plot(M.Time_s_(iT:nT:nt*nT), p_mean(iT:nT:end),'Color',cols(3,:),'LineWidth',2)
    plot(powSOWFA(iT:nT:end,2), powSOWFA(iT:nT:end,4),'-k','LineWidth',1)
    hold off
    ylim([minP-.1, maxP + .1])
    xlim([M.Time_s_(1), M.Time_s_(end)])
    grid on
    title(['Power generated, T' num2str(iT-1)])
    ylabel('Power in MW')
    xlabel('Time in s')
end

% Predicted vs actual power


end

% function private_power(path)
% data_dirs = dir([path filesep 'EnKF_Power']);
% 
% end