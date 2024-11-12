addpath(genpath("/Users/marcusbecker/surfdrive/PhD_Surf/02_Communication/00_Plot_lib"))

%% Create a test environment

% Number of turbines
T.nT = 10;

% x,y,z location of the turbines
T.posBase = [
    2681.25712384114,581.073357148974,0
    2275.94463924756,1269.94121808646,0
    3170.98183442861,1341.08248067317,0
    1870.63215465397,1958.80907902393,0
    3660.55143066563,2102.35490920341,0
    1465.31967006039,2647.67693996141,0
    3939.99281453319,3240.48494262761,0
    1060.00718546681,3336.54480089889,0
    2997.06559432448,3823.37541895182,0
    2136.36407247608,4418.92664285103,0];

T.nT = size(T.posBase,1);

% Wind (direction) state (minial working example)
T.States_WF      = zeros(T.nT,3);
T.States_WF(:,2) = 270;
T.StartI         = 1:T.nT;

% Turbine diameter
T.D = ones(T.nT,1)*178.4;

% Free wind speed
Uinfty  = 8;

% Action horizon length simutaneously the number of optimization time steps
ah      = 20;
dof     = 2;

% Time step in s
deltaT  = 5;

%% RUN THE ALGORITHM
[T_star,combined,ah,ph,solve_order,t2o, t2g, opt_ph] = ...
    calcT2(T, Uinfty, ah, dof, deltaT);

%% DEMO Optimization execution
% Optimization variables for all turbines
%   The assumption is that x is a (ah*nT x 1) vector of optimization
%   variables that need to be set. (ah) are the action horizon steps, nT
%   the number of turbines.
x = zeros(dof*T.nT,1);
t = 0:T.nT-1;

% ///// Plotting
x_history = zeros(length(x),size(t2o,1));
x_hist_i  = 1;
T_star2t = repmat(t,ph,1); T_star2t = T_star2t(:)';

% Solve each independent group
for iG = 1:length(solve_order)      % Can be executed in parallel threads

% Solve dependent optimization problems within the group
for iO = 1:length(solve_order{iG})  % Can NOT be executed in parallel
    subset_x    = repmat(combined(solve_order{iG}(iO),:)>0,dof,1);
    subset_t    = t2g(iG,:)>0;

    % Get subset of x that is optimized in this simulation
    x_opt = x(subset_x(:));

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
        'actuating turbine(s) ' num2str(unique(t(combined(solve_order{iG}(iO),:)>0))')...
        ', optimizing ' num2str(length(x_opt)*dof) ' values across ' ...
        num2str(opt_ph(solve_order{iG}(iO))) ' time steps'])
    
    % Cost function setup
    %   Calculate the power generated of the reduced wind farm based on 
    %   x_opt over the reduced prediction horizon  (DUMMY VALUE)
    p = rand(length(t_opt)*ph_opt,1); 
    % Combine the power generated as a cost function
    J = -T_star_opt * p;
    % Get an optimal x (DUMMY VALUE)
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
xlabel('Optimization Nr.')
ylabel('x vector elements')
colormap(WiBk(length(solve_order)+1))
grid on

clear x_history iG x2t t_opt x_hist_i x_opt subset_x subset_t x t iO
