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

function time_matrix = calcTimes(T, Uinfty)
%CALCTIMES calculates the time between the turbines based on the free wind
%speed and direction

% http://dx.doi.org/10.1016/j.renene.2017.06.065
U_adv = 0.7*Uinfty;

% Time matrix lists the time it takes for the actuation signal at T_col to
% reach T_row
time_matrix = zeros(T.nT, T.nT);

for iT = 1:T.nT
    dist    = T.posBase(:,1:2)' - T.posBase(iT,1:2)';
    phi     = angSOWFA2world(T.States_WF(T.StartI(iT),2));
    dist    = [cos(phi), sin(phi); -sin(phi), cos(phi)] * dist;

    time_matrix(iT,:) = dist(1,:)./U_adv;
    
    % Remove entries that are upstream (-> -1), too far away cross stream
    % (-> -2)
    time_matrix(iT,time_matrix(iT,:)<0) = -1;
    time_matrix(iT,and(time_matrix(iT,:)>0, abs(dist(2,:)./T.D(iT))>2)) = -2;
end


end

