function [tr, g0] = mr_2_tr(tn, x, cu, rn, g0)
%MR_1_TR Max-rate 1-dof yaw trajectory
arguments
    tn (:,1) double % Time series, has to be equal to/between 0 and 1
    x  (1,:) double % Optimization var, has to be equal to/between 0 and 1
    cu (1,1) double % Ratio of control update to time window of length 1
    rn (1,1) double % Normalised rate limit (to time window of length 1)
    g0 (1,:) double = 0 % Current angle
end

x = reshape(x,2,[]);
% Saturation function
sat = @(t) min(1,max(0,t));

% Duration of the change based on the output gamma
t_d = @(g_so) abs((g_so-.5));

% Starting time
t_s = @(t_so,g_so) t_so.*(1-t_d(g_so));

% Calculate trajectory
tr = 2*(x(1,:)-.5) .* rn .* ...
    sat((tn-t_s(x(2,:),x(1,:)))./(2*t_d(x(1,:)))) + g0;

if cu>1
    % Update outside of the trajectory
    % Calculate the offset at the control update
    g0 = sign(x(1,:)-.5).*rn.*t_d(x(1,:)).*sat((1)./t_d(x(1,:))) + g0;

elseif cu>0
    % Update inside of the trajectory
    % Calculate the offset at the control update
    g0 = sign(x(1,:)-.5).*rn.*t_d(x(1,:)).*...
        sat((cu-t_s(x(2,:),x(1,:)))./t_d(x(1,:))) + g0;
else
    error('The ratio of control update to time window is negative')
end
end

