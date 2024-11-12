addpath(genpath("/Users/marcusbecker/surfdrive/PhD_Surf/02_Communication/00_Plot_lib"))

%% Create a test environment

% Number of turbines
T.nT = 10;

% x,y,z location of the turbines
T.posBase = [
     804.179581909827,  1524.4079389476,    0
    1374.53115876326,   2084.33312929074,   0
    1660.08831670344,   1233.09299304884,   0
    1944.88273561669,   2644.25831963388,   0
    2517.1853053877,     942.23417508319,   0
    2515.23431247012,   3204.18350997702,   0
    3689.1110020598,    946.431969696743,   0
    3085.58588932356,   3764.10870032017,   0
    4026.57239959874,   2002.36419061391,   0
    4396.2106632155,    2981.57607987506,   0];

% Wind (direction) state (minial working example)
T.States_WF      = zeros(T.nT,3);
T.States_WF(:,2) = 225;
T.StartI         = 1:T.nT;

% Turbine diameter
T.D = ones(T.nT,1)*178.4;

% Free wind speed
Uinfty  = 8;

% Action horizon length simutaneously the number of optimization time steps
ah      = 5;

% Time step in s
deltaT  = 5;

%% RUN THE ALGORITHM
[T_star,combined,ah,ph,solve_order,t2o, t2g, opt_ph] = calcT(T, Uinfty, ah, deltaT);

%% DEMO Optimization execution
% Optimization variables for all turbines
%   The assumption is that x is a (ah*nT x 1) vector of optimization
%   variables that need to be set. (ah) are the action horizon steps, nT
%   the number of turbines.
x = zeros(ah*T.nT,1);
t = 0:T.nT-1;

% ///// Plotting
x_history = zeros(length(x),size(t2o,1));
x_hist_i  = 1;
x2t = repmat(t,ah,1); x2t = x2t(:);
T_star2t = repmat(t,ph,1); T_star2t = T_star2t(:)';

% Solve each independent group
for iG = 1:length(solve_order)      % Can be executed in parallel threads

% Solve dependent optimization problems within the group
for iO = 1:length(solve_order{iG})  % Can NOT be executed in parallel
    subset_x    = combined(solve_order{iG}(iO),:)>0;
    subset_t    = t2g(iG,:)>0;

    % Get subset of x that is optimized in this simulation
    x_opt = x(subset_x);

    % Get subset of the turbines that are simulated in this simulation
    t_opt = t(subset_t);

    % Get the prediction horizon of the optimization
    ph_opt = opt_ph(solve_order{iG}(iO));

    % Get the vector the reduced simulation output needs to be multiplied
    % by:
    T_star_opt = reshape(T_star(solve_order{iG}(iO),:),ph,[]);
    %   Dim 1: Shorten to prediction horizon
    %   Dim 2: Select only relevant turbines
    T_star_opt = T_star_opt(1:ph_opt, t2g(iG,:)>0);
    T_star_opt = T_star_opt(:)';

    
    % ============== Optimize ============== %
    disp(['Simulating the turbine(s) ' num2str(t_opt) ', ' ...
        'actuating turbine(s) ' num2str(unique(x2t(combined(solve_order{iG}(iO),:)>0))')...
        ', optimizing ' num2str(length(x_opt)) ' values across ' ...
        num2str(opt_ph(solve_order{iG}(iO))) ' time steps'])
    
    % Cost function setup
    %   Calculate the power generated of the reduced wind farm based on 
    %   x_opt over the reduced prediction horizon 
    p = rand(length(t_opt)*ph_opt,1); 
    % Combine the power generated as a cost function
    J = -T_star_opt * p;
    % Get an optimal x
    x_opt(:) = iG;
    % ====================================== %


    % Store results in optimization variable vector
    x(subset_x) = x_opt;

    
    % ///// Plotting
    x_history(:,x_hist_i) = x;
    x_hist_i = x_hist_i+1;
end
end

% ///// Plotting
figure
imagesc(x_history)
xlabel('Optimizations')
ylabel('x vector elements')
colormap(WiBk(length(solve_order)+1))
grid on

clear x_history iG x2t t_opt x_hist_i x_opt subset_x subset_t x t iO

%% Overview plot of the relation between the optimizations, the turbines 
% and runtime
figure
subplot(1,3,1)
imagesc(0:T.nT-1,[],t2o)
title("Actuated turbines in optimizations")
colormap(WiBk(2))
xlabel('Turbines')
ylabel('Optimization')

subplot(1,3,2)
y_data = repmat(.5:1:length(opt_ph)+.5,2,1);
x_data = repmat(opt_ph',2,1);
patch([0;x_data(:);0], y_data(:),'k')
hold on
plot([1,1]*ah, [.5, length(opt_ph)+.5], '--r', 'LineWidth',2)
xlabel('Time in (s)')
ylabel('Optimization')
set(gca,'YDir','reverse')
ylim([.5, length(opt_ph)+.5])
title('Simulated time')

subplot(1,3,3)
tsim = zeros(size(t2o));
for iG = 1:length(solve_order)
    for iO = 1:length(solve_order{iG})
        tsim(solve_order{iG}(iO),:) = t2g(iG,:);
    end
end
imagesc(0:T.nT-1,[],tsim)
title("Simulated turbines in optimizations")
colormap(WiBk(2))
xlabel('Turbines')
ylabel('Optimization')

%%
figure
imagesc(0:T.nT-1,[],t2g)
title("Turbines in group (t2g)")
colormap(Bu(2))

%%
figure
imagesc(T_star)
title("T^*")
colormap(Bu(2))

%% Plot turbine group solve order
i_map = 0;

figure
hold on
for iT = 1:T.nT
    plot([0 10000*cos(angSOWFA2world(T.States_WF(iT,2)))] + T.posBase(iT,1), ...
        [0 10000*sin(angSOWFA2world(T.States_WF(iT,2)))] + T.posBase(iT,2), ...
        '-k','LineWidth',1)
    
end


for iG = 1:length(solve_order)
    if mod(i_map,3) == 0
        cols = flipud(Bu(length(solve_order{iG})+2));
    elseif mod(i_map,3) == 1
        cols = (Rd(length(solve_order{iG})+2));
    else
        cols = WiBk(length(solve_order{iG})+2);
    end
    cols(1,:) = [];
    cols(end,:) = [];
    i_map = i_map+1;
    for iO = 1:length(solve_order{iG})
        indO    = solve_order{iG}(iO);
        sz      = 100 * (length(solve_order{iG}) - iO + 1);
        x_T     = T.posBase( t2o(indO,:)>0 , 1);
        y_T     = T.posBase( t2o(indO,:)>0 , 2);
        scatter(x_T, y_T, sz, cols(iO,:),"filled")
    end
end

for iT = 1:T.nT
    text(T.posBase(iT,1),T.posBase(iT,2),['T' num2str(iT-1)], ...
        'FontWeight', 'bold')
    
end


hold off
axis equal
xlim([0,5000])
ylim([0,5000])

clear sz i_map cols iG iO x_T y_T





%% Run demo of the optimization
% Creates a plot visualization of what's happening
%
% Optimization variables for all turbines
%   The assumption is that x is a (ah*nT x 1) vector of optimization
%   variables that need to be set. (ah) are the action horizon steps, nT
%   the number of turbines.
x = zeros(ah*T.nT,1);
t = 0:T.nT-1;

% ///// Plotting
x_history = zeros(length(x),size(t2o,1));
x_hist_i  = 1;
x2t = repmat(t,ah,1); x2t = x2t(:);
T_star2t = repmat(t,ph,1); T_star2t = T_star2t(:)';

phi = deg2rad(270 - mean(T.States_WF(:,2)));
cols = flipud(RdBk(8));
% Solve each independent group
for iG = 1:length(solve_order)      % Can be executed in parallel threads

% Solve dependent optimization problems within the group
for iO = 1:length(solve_order{iG})  % Can NOT be executed in parallel
    subset_x    = combined(solve_order{iG}(iO),:)>0;
    subset_t    = t2g(iG,:)>0;

    % Get subset of x that is optimized in this simulation
    x_opt = x(subset_x);

    % Get subset of the turbines that are simulated in this simulation
    t_opt = t(subset_t);

    % Get the prediction horizon of the optimization
    ph_opt = opt_ph(solve_order{iG}(iO));

    % Get the vector the reduced simulation output needs to be multiplied
    % by:
    T_star_opt = reshape(T_star(solve_order{iG}(iO),:),ph,[]);
    %   Dim 1: Shorten to prediction horizon
    %   Dim 2: Select only relevant turbines
    T_star_opt = T_star_opt(1:ph_opt, t2g(iG,:)>0);
    T_star_opt = T_star_opt(:)';

    T_measure = reshape(T_star_opt, [], length(t_opt));


    Tmp_posBase = T.posBase(subset_t,:);

    %points = zeros(ength(x_opt), 2);
    combined_tmp = reshape(combined(solve_order{iG}(iO),:), ah, []);
    combined_tmp = combined_tmp(:, t_opt+1);

    points_x = nan(size(combined_tmp));
    points_y = nan(size(combined_tmp));
    
    figure
    for i_ph = 1:opt_ph(solve_order{iG}(iO))
        for iT = 1:size(Tmp_posBase,1)
            plot([0, cos(phi)]*10000  + Tmp_posBase(iT,1), ...
                [0, sin(phi)]*10000 + Tmp_posBase(iT,2),'--k')
            if iT == 1; hold on; end
        end
        
        scatter(T.posBase(:,1), T.posBase(:,2),20,[.3 .3 .3], ...
            "filled", "x", 'MarkerEdgeColor','flat')

        % Set starting point to turbine
        if i_ph <= ah
            points_x(i_ph,:) = Tmp_posBase(:,1);
            points_y(i_ph,:) = Tmp_posBase(:,2);
        end


        points_x = points_x + cos(phi)*Uinfty*0.7396*deltaT;
        points_y = points_y + sin(phi)*Uinfty*0.7396*deltaT;

        scatter(points_x(combined_tmp<1), points_y(combined_tmp<1), ...
            50,cols(2,:),"filled")

        scatter(points_x(combined_tmp>0), points_y(combined_tmp>0), ...
            50,cols(end-1,:),"filled")
        
        
        % Plot involved turbines
        

        scatter(Tmp_posBase( T_measure(i_ph,:)>0,1), ...
            Tmp_posBase(T_measure(i_ph,:)>0,2), 80, ...
              'k','filled', 'o', 'MarkerEdgeColor', 'flat' )

        scatter(Tmp_posBase( T_measure(i_ph,:)==0,1), ...
            Tmp_posBase(T_measure(i_ph,:)==0,2), 60, ...
              'w','filled', 'o', 'MarkerEdgeColor', 'k' )

        hold off
        axis equal
        xlim([0,5000])
        ylim([0,5000])
        grid on
        title(['Actuated turbine(s) ' num2str(unique(x2t(combined(solve_order{iG}(iO),:)>0))') ', Timestep ' num2str(i_ph)])
        colormap(cols(2:end-1,:))
        pause(.5)
    end
    

    % ============== Optimize ============== %
    disp(['Simulating the turbine(s) ' num2str(t_opt) ', ' ...
        'actuating turbine(s) ' num2str(unique(x2t(combined(solve_order{iG}(iO),:)>0))')...
        ', optimizing ' num2str(length(x_opt)) ' values across ' ...
        num2str(opt_ph(solve_order{iG}(iO))) ' time steps'])
    
    % Cost function setup
    %   Calculate the power generated of the reduced wind farm based on 
    %   x_opt over the reduced prediction horizon 
    p = rand(length(t_opt)*ph_opt,1); 
    % Combine the power generated as a cost function
    J = -T_star_opt * p;
    % Get an optimal x
    x_opt(:) = iG;
    % ====================================== %


    % Store results in optimization variable vector
    x(subset_x) = x_opt;

    
    % ///// Plotting
    x_history(:,x_hist_i) = x;
    x_hist_i = x_hist_i+1;
end
end

% ///// Plotting
figure
imagesc(x_history)
xlabel('Optimizations')
ylabel('x vector elements')
colormap(WiBk(length(solve_order)+1))
grid on

clear x_history iG x2t t_opt x_hist_i x_opt subset_x subset_t x t iO

