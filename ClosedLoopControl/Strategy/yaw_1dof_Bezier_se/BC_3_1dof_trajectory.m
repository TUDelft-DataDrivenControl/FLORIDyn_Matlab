% Copyright (C) <2024>, M Becker
%
% List of the contributors to the development of FLORIDyn: see LICENSE file.
% Description and complete License: see LICENSE file.
	
% This program (FLORIDyn) is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program (see COPYING file).  If not, see <https://www.gnu.org/licenses/>.
% ======================================================================= %
% Updated: 16. Dez. 2024, M. Becker
% ======================================================================= %

function [tr, g0, gd] = BC_3_1dof_trajectory(t_n,x,gd,rn,cu,g0)
%BC_3_1DOF_TRAJECTORY Derives the new trajectory based on the optimization
%value and the current derivative. The function returns the trajectory and
%the offset & derivative at the time of the next control update.

arguments
    t_n (:,1) double % Time series, has to be equal to/between 0 and 1
    x   (1,:) double % Optimization var, has to be equal to/between 0 and 1
    gd  (:,1) double % Normalised derivative (to time window of length 1)
    rn  (1,1) double % Normalised rate limit (to time window of length 1)
    cu  (1,1) double % Ratio of control update to time window of length 1
    g0  (1,:) double = 0 % Current angle
end

% Convert x into the gamma 2 point
g2 = bc_3_o1_to_g2(x, gd/3, rn);

% Calculate the trajectory
tr = bc_3_tr(t_n,0,gd/3,g2,g2) + g0;

if cu>1
    % Update outside of the trajectory
    % Calculate the offset at the control update
    g0 = bc_3_tr(1,0,gd/3,g2,g2) + g0;
    
    % Calculate the derivative at the control update
    gd = bc_3_d(1,0,gd/3,g2,g2); % == 0
elseif cu>0
    % Update inside of the trajectory
    % Calculate the offset at the control update
    g0 = bc_3_tr(cu,0,gd/3,g2,g2) + g0;
    
    % Calculate the derivative at the control update
    gd = bc_3_d(cu,0,gd/3,g2,g2);
else
    error('The ratio of control update to time window is negative')
end


end

