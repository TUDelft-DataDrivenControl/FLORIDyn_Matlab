clear all
addpath('./WindField/Velocity_Interpolation_wErrorCov/');

% Simple linear data: speed = 5 at t=0, speed = 15 at t=10
times = [0.0, 10.0];
speeds = [5.0, 15.0];
data = [times(:), speeds(:)];  % Column-wise concatenation

% Simple 2x2 identity Cholesky: noise still generated but uncorrelated
chol = chol(eye(2));  % chol returns upper triangular by default in MATLAB

% Define a WindVelMatrix function equivalent
wind = WindVelMatrix(data, chol);

% Test 5: Single turbine (scalar iT input)
iT_scalar = 1;
vel_scalar = getWindSpeedT(wind, iT_scalar, 5.0);

% assert(length(vel_scalar) == 1);
% assert(abs(vel_scalar(1) - 10.0) <= 3.0);


%% --- Supporting functions ---

function wind = WindVelMatrix(data, chol)
    % Store data and chol in a struct (MATLAB equivalent of a Julia type)
    wind.Data = data;
    wind.CholSig = chol;
end

% function vel = getWindSpeedT(wind, iT, tQuery)
%     % Simple linear interpolation to get speed at time tQuery for turbine iT
%     % Assuming wind.data(:,1) = times, wind.data(:,2) = speeds
% 
%     times = wind.data(:,1);
%     speeds = wind.data(:,2);
% 
%     % Interpolate speed linearly at time tQuery
%     speed_t = interp1(times, speeds, tQuery, 'linear');
% 
%     % Add noise using chol matrix (noise generation simplified)
%     noise = wind.chol * randn(2,1);
% 
%     % For this example just take the first element noise as perturbation
%     speed_noisy = speed_t + noise(1);
% 
%     % Return speed as scalar (for single turbine)
%     vel = speed_noisy;
% end
