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

function YawDataNew = condenseSOWFAYaw(YawData)
%CONDENSESOWFAYAW Reduce yaw data to what is needed
diff = sum(...
    abs(YawData(2:end-1,2:end) - YawData(1:end-2,2:end)) + ...
    abs(YawData(2:end-1,2:end) - YawData(3:end,2:end)),2);
ind_important = [1;find(diff>0)+1;size(YawData,1)];
YawDataNew = YawData(ind_important,:);
end

