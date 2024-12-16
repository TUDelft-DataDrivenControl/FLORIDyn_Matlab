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

function T = correctTi(T,Wind,SimTime)
%CORRECTTI Correction of the turbulent intensity
%% Get Data
TI = getDataTI(Wind,T,SimTime);
%% Correct
for iT = 1:T.nT
    if isempty(T.dep{iT})
        T.States_WF(T.StartI(iT),3) = TI(iT);
        continue
    end
    
    % Assign amb. TI based on what the influencing OPs carry
    if size(T.intOPs{iT},1)==1
        % Only one turbine influencing
        T.States_WF(T.StartI(iT),3) = ...
            T.States_WF(T.intOPs{iT}(1),3) * T.intOPs{iT}(2) + ...
            T.States_WF(T.intOPs{iT}(3),3) * T.intOPs{iT}(4);
    else
        % More than one, currently only taking the mean, but could be more
        % sophisticated (e.g. only from uninfluenced wind turbines,
        % weighted by distance to OP)
        T.States_WF(T.StartI(iT),1) = mean(...
            T.States_WF(T.intOPs{iT}(:,1),3) .* T.intOPs{iT}(:,2) + ...
            T.States_WF(T.intOPs{iT}(:,3),3) .* T.intOPs{iT}(:,4));
    end
end
end

