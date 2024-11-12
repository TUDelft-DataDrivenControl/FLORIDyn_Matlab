function [tr, g0] = mr_1_tr(tn, x, cu, rn, g0)
%MR_1_TR Max-rate 1-dof yaw trajectory
arguments
    tn (:,1) double % Time series, has to be equal to/between 0 and 1
    x  (1,:) double % Optimization var, has to be equal to/between 0 and 1
    cu (1,1) double % Ratio of control update to time window of length 1
    rn (1,1) double % Normalised rate limit (to time window of length 1)
    g0 (1,:) double = 0 % Current angle
end

% Saturation function
sat = @(t) min(1,max(0,t));

% Duration of the change based on the output gamma
t_d = @(g_so) abs((g_so-.5));

% Calculate trajectory
tr = sign(x-.5).*rn.*t_d(x).*sat((tn-0)./t_d(x)) + g0;

if cu>1
    % Update outside of the trajectory
    % Calculate the offset at the control update
    g0 = sign(x-.5).*rn.*t_d(x).*sat((1)./t_d(x)) + g0;

elseif cu>0
    % Update inside of the trajectory
    % Calculate the offset at the control update
    g0 = sign(x-.5).*rn.*t_d(x).*sat((cu-0)./t_d(x)) + g0;
else
    error('The ratio of control update to time window is negative')
end
end

