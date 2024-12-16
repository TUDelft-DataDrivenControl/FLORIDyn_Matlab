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

function storeCLInfo(T,Sim,it)
% STORECLINFO writes the current OP position in the down / crosssteam
% coordinate system of their respective turbine.

for iT =1:T.nT
    % Calculate dw, cw of OP chain & T pos in K1 of iT
    if iT<T.nT
        delta_x_op = ...
            T.States_OP(T.StartI(iT):T.StartI(iT+1)-1, 1) - ...
            T.posBase(iT,1);
        delta_y_op = ...
            T.States_OP(T.StartI(iT):T.StartI(iT+1)-1, 2) - ...
            T.posBase(iT,2);
    else
        delta_x_op = ...
            T.States_OP(T.StartI(iT):end, 1) - ...
            T.posBase(iT,1);
        delta_y_op = ...
            T.States_OP(T.StartI(iT):end, 2) - ...
            T.posBase(iT,2);
    end

    % Add other turbine loactions as positions 
    for iiT = 1:T.nT
        %if iiT == iT; continue; end
        delta_x_op(end+1) = ...
            T.posBase(iiT,1) - T.posBase(iT,1);
        delta_y_op(end+1) = ...
            T.posBase(iiT,2) - T.posBase(iT,2);
    end
    phi = angSOWFA2world(T.States_WF(T.StartI(iT),2));

    % Rotate world positions into down-/crosssteam positions of turbine iT
    delta_dw_op = cos(phi)*delta_x_op + sin(phi)*delta_y_op;
    delta_cw_op = - sin(phi)*delta_x_op + cos(phi)*delta_y_op;

    % Write data (create new file in first iteration)
    if it == 1
        writematrix(delta_dw_op, ...
            [Sim.PathToSim filesep 'Results' filesep 'dw_to_T' ...
            num2str(iT-1) '.csv'],...
            'WriteMode','overwrite')
        writematrix(delta_cw_op, ...
            [Sim.PathToSim filesep 'Results' filesep 'cw_to_T' ...
            num2str(iT-1) '.csv'],...
            'WriteMode','overwrite')
    else
        writematrix(delta_dw_op, ...
            [Sim.PathToSim filesep 'Results' filesep 'dw_to_T' ...
            num2str(iT-1) '.csv'],...
            'WriteMode','append')
        writematrix(delta_cw_op, ...
            [Sim.PathToSim filesep 'Results' filesep 'cw_to_T' ...
            num2str(iT-1) '.csv'],...
            'WriteMode','append')
    end
end