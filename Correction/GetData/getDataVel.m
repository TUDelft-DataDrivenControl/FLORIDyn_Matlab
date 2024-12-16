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

function [U, Wind] = getDataVel(Wind,T,SimTime,tmpM,paramFLORIS)
%GETDATAVEL retreves the data for the wind speed for every turbine

switch Wind.Input_Vel
    case 'I_and_I'
        [U,Wind.Vel] = getWindSpeedT(Wind.Vel,(1:T.nT)',SimTime,...
            T.States_WF(T.StartI,2),paramFLORIS.p_p);
        if strcmp(Wind.Input_Vel,'I_and_I')
            if (SimTime-Wind.Vel.StartTime)>Wind.Vel.WSE.Offset
                % Ufree = Ueff/reduction
                U = U./tmpM(:,1);
            end
        end
    case 'ZOH_wErrorCov'
        U = getWindSpeedT(T.States_WF(T.StartI,1), Wind.Vel.ColSig);
    case 'RW_with_Mean'
        U = getWindSpeedT(T.States_WF(T.StartI,1),Wind.Vel);
    case 'EnKF_InterpTurbine'
        U = getWindSpeedT_EnKF(Wind.Vel,(1:T.nT)',SimTime);
    case 'EnKF_RW'
        U = T.States_WF(T.StartI,1) ;
        U = U + (randn(1,length(U))*Wind.Vel.CholSig)';
    case 'EnKF_ZOH'
        U = T.States_WF(T.StartI,1);
    case 'CLC_weighted_ZOH'
        % Apply weighted states on the turbine
        U = T.C_Vel*T.States_WF(:,1);
    otherwise
        U = getWindSpeedT(Wind.Vel,(1:T.nT)',SimTime);
end

end

