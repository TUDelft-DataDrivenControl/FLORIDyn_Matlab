function [t,yaw] = BC_3_2dof_trajectory(t, t_res, yaw_n, yaw_0)
%BC_3_2DOF_TRAJECTORY Converts yaw setpoints into a trajectory
%
% INPUTS
%   t_n    [nT x 1] Time of the normalised trajectory
%   yaw_n  [nT x 2] Yaw setpoints
%   yaw_0  [nT x 1] Yaw orientation
%
% OUTPUTS
%   t      [nt x 1]  Time stamps for the steps of all turbines
%   yaw    [nt x nT] Yaw orientation of the turbines
%%
nT = length(t);
nt = round(max(t)/t_res)+1;
t_n = linspace(0,1,nt)';
yaw = zeros(nt,nT);

% should be doable without 'for', but I have no time right now
for iT = 1:nT
    yaw(:,iT) = 3*(1-t_n).^2.*t_n*yaw_n(iT,1) + ...
        yaw_n(iT,2)*(3*(1-t_n).*t_n.^2 + t_n.^3) + yaw_0(iT);

end

t = t_n * max(t);

end

