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


function cov = GaspariAndCohn1999(l,d)
%GASPARIANDCOHN1999
% "To achieve covariance localization by Schur product a function, ρ, 
% is normally defined to be a correlation function with local support. 
% Local support is a term meaning that the function is only non zero in a 
% small (local) region and is zero elsewhere. The correlation function is 
% commonly taken to be the compactly supported 5th order piecewise rational
% function as defined in Gaspari and Cohn (1999)"
%       Localization in the ensemble Kalman Filter, R.E. Petrie, 2008
% ======================================================================= %
% l = cut off length
% d = distance(s) to the measurement location
% ======================================================================= %
c = sqrt(10/3)*l;
    
cov = .5*(sign(d) - sign(d-c)).*(...
        (((-1/4*(d./c) + 1/2).*(d./c) + 5/8).*(d./c) - 5/3) ...
        .*(d./c).*(d./c) + 1)...
        + .5*(sign(d-c) - sign(d-2*c)).*(...
        ((((1/12*(d./c) - 0.5).*(d./c) + 5/8).*(d./c) + 5/3).*(d./c) - 5).*(d./c) ...
        + 4 - 2/3*(c./d));
cov(isnan(cov))=1;
end

%% Test plot:
% figure
% x = linspace(0,50);
% l = 10;
% plot(x,GaspariAndCohn1999(10,x))
% hold on
% scatter(l,GaspariAndCohn1999(l,l),'filled')
% hold off
% grid on
