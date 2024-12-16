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

function [d, EnKF] = EnKF_GetWFInputs(EnKF,Sim,T,paramFLORIS)
%ENKF_GETWFINPUTS reads inputs from a validation source and returns them as
%a measurement vector

for iW = 1:EnKF.nS+1
    SimTime = Sim.StartTime+(iW-1)*Sim.TimeStep;
    
    %% Wind speed
    if strcmp(EnKF.Wind.Input_Vel,'I_and_I')
        [U,EnKF.Wind.Vel] = getWindSpeedT(EnKF.Wind.Vel,(1:T.nT)',SimTime,...
            T.States_WF(T.StartI,2),paramFLORIS.p_p);
        if (SimTime-EnKF.Wind.Vel.StartTime)>EnKF.Wind.Vel.WSE.Offset
            % Ufree = Ueff/reduction
            %U = U./tmpM(:,1);
%             warning(['The EnKF can not yet use the I&I as intended, '...
%                 'proceeds with U_eff=U_free'])
        end
    elseif strcmp(EnKF.Wind.Input_Vel,'ZOH_wErrorCov')
        U = getWindSpeedT(U, EnKF.Wind.Vel.ColSig);
    elseif strcmp(EnKF.Wind.Input_Vel,'RW_with_Mean')
        U = getWindSpeedT(U,EnKF.Wind.Vel);
    else
        U = getWindSpeedT(EnKF.Wind.Vel,(1:T.nT)',SimTime);
    end
    
    %% Wind direction
    if strcmp(EnKF.Wind.Input_Vel,'RW_with_Mean')
        phi = getWindDirT(phi,EnKF.Wind.Dir);
    else
        phi = getWindDirT(EnKF.Wind.Dir,(1:T.nT)',SimTime);
    end
    
    %% Ambient turbulence intensity
    TI = getWindTiT(EnKF.Wind.TI,1:T.nT,...
        SimTime)';
end
d = [U,phi,TI];
end

