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

%% Before plotting
figure
if dark_mode; set(gcf,'Color','k'); end

%% Isosurface plotting example
s = isosurface(X,Y,Z,U,8);
p = patch(s);
set(p,'FaceColor',[1,1,1]);
set(p,'FaceAlpha',0.1);

%% After plotting
ax = gca;
ax.YColor = 'w';
ax.XColor = 'w';
ax.Color = 'k';
ax.GridColor = 'w';