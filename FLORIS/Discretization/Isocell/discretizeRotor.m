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

function [RPs,w] = discretizeRotor(nRP)
%DISCRETIZEROTOR discretizes the rotor plane into nRP segments. The
% algorithm returns the normalized center location \in [-.5,.5] and the
% relative area the segment represents.

%% The isocell algorithm 
% It faces certain limitations and aims to achieve nRP cells but might end
% up with a slightly different number. For details, see the publication of
% Masset et al.:
%  https://orbi.uliege.be/bitstream/2268/91953/1/masset_isocell_orbi.pdf
% We choose N1 = 3 here, 4 or 5 are also viable options, 3 is close to
% optimal
N1 = 3;
n = round(sqrt(nRP/N1));
%eRP = nRP - N1*n^2; % Difference from calculated to actually applied


%% Generate pattern
% Radial thickness of each ring
dltR  = 1/n;
nC      = N1*n^2;
% Each ring i has (2i-1)*N1 cells
RPs = zeros(nC,3);
for i = 1:n
    % Segments in ring
    nR = (2*i-1) * N1;
    
    % Start index and end index
    i_e = sum((2*(1:i)-1)*N1);
    i_s = i_e - nR + 1;
    
    phi = (1:nR)/nR*2*pi;
    RPs(i_s:i_e,2) = 0.5*cos(phi)*dltR*(.5+(i-1));
    RPs(i_s:i_e,3) = 0.5*sin(phi)*dltR*(.5+(i-1));
end

w = ones(nC,1)*1/nC;

% figure
% scatter3(RPs(:,1),RPs(:,2),RPs(:,3))
% axis equal
end

