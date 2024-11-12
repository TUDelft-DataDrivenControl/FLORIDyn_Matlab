function [tr, g0, gd] = BC_3_2dof_trajectory(t_n,x,gd,rn,cu,wa,g0)
%BC_3_1DOF_TRAJECTORY Derives the new trajectory based on the optimization
%value and the current derivative. The function returns the trajectory and
%the offset & derivative at the time of the next control update.

arguments
    t_n (:,1) double % Time series, has to be equal to/between 0 and 1
    x   (:,1) double % Optimization var, has to be equal to/between 0 and 1
    gd  (:,1) double % Normalised derivative (to time window of length 1)
    rn  (1,1) double % Normalised rate limit (to time window of length 1)
    cu  (1,1) double % Ratio of control update to time window of length 1
    wa  (:,1) double % Wind angle at the end of the time window (normed)
    g0  (1,:) double = 0 % Current angle
end

x = reshape(x,[],2);

% Convert x into the gamma 2 point
g2 = x(:,1)*rn/3 - rn/6 + gd/3;

% Convert gamma 2 & wind angle into g3
g3 = 2*x(:,2)*rn/3 - rn/3 + g2;

% Calculate the trajectory
tr = bc_3_tr(t_n,0,gd/3,g2,g3) + g0;

if cu>1
    % Update outside of the trajectory
    % Calculate the offset at the control update
    g0 = bc_3_tr(1,0,gd/3,g2,g3) + g0;
    
    % Calculate the derivative at the control update
    gd = bc_3_d(1,0,gd/3,g2,g3); % == 0
elseif cu>0
    % Update inside of the trajectory
    % Calculate the offset at the control update
    g0 = bc_3_tr(cu,0,gd/3,g2,g3) + g0;
    
    % Calculate the derivative at the control update
    gd = bc_3_d(cu,0,gd/3,g2,g3);
else
    error('The ratio of control update to time window is negative')
end


end

