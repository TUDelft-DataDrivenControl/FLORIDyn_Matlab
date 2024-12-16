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

% CLC EnKF script

% Set inputs
Sim.StartTime   = time.current - time.stepEnKF;
Sim.EndTime     = time.current;
Sim.nSimSteps   = EnKF.nS;

for iE = 1:EnKF.nE
    % ===== Assign relevant ensemble states
    T = EnKF_AssignEnStates(EnKF,T,iE);
    
    % =========== Run simulation ================
    [T,M,Vis] = ...
        FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
    
    % ===== Stores relevant measurements and states
    EnKF = EnKF_StoreEnStates(EnKF,M, T, iE);
end



%% Combination
% Calculate "true" state & state error covariance
if EnKF.Vel.Correct
    % ===== Projection on the mean "true" state
    %           Velocity
    [EnKF,~]        = EnKF_projectOntoTrueState(EnKF,Sim,T,1);
    %           Direction
    [EnKF, truePos] = EnKF_projectOntoTrueState(EnKF,Sim,T,2);

    % Calc C_phi and C_u based on true states
    [C_u, C_phi] = EnKF_calcC_uC_phi(truePos,...
        mean(EnKF.States.Dir,2), Sim.Dyn, T.posBase, T.nOP,...
        Sim.TimeStep);
    
    % ===== State Error covariance matrix
    C_xx_Vel = EnKF_CalcStateCov(EnKF.States.Vel,EnKF.nE,T.StartI);

    % ===== Power based on projected states
    P = EnKF_calcPower(C_u, EnKF, T.D, paramFLORIS, T.StartI);
    
    % ===== State-to-Output & Output error covariance matrices
    [C_xy_Vel, C_yy_Vel] = ...
        EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE,...
         'Method', 'InputOutput', 'Output',P);%EnKF.Output.Pow);
end

if EnKF.Dir.Correct
    C_xx_Dir = EnKF_CalcStateCov(EnKF.States.Dir,EnKF.nE,T.StartI);
end

if EnKF.TI.Correct
    C_xx_TI = EnKF_CalcStateCov(EnKF.States.TI,EnKF.nE,T.StartI);
end

%% Get measurements from validation source
%   Calls the functions to get wind speed, direction and amb.
%   turbulence as the normal simulation would and gets [U,phi,TI]
[d, EnKF] = EnKF_GetWFInputs(EnKF,Sim,T,paramFLORIS);

% Calculate reduction of free wind speed at the turbines & apply to meas.
% r = EnKF_GetAveForeignReduction(EnKF.M,T.nT); 
% d(:,1) = d(:,1)./(r*10^(-2)); 

d_P = interp1(powSOWFA(1:T.nT:end,2),reshape(powSOWFA(:,3),T.nT,[])',Sim.EndTime)';


%% Correct Ensemble states
% [1] Eq.4.37
% x_a,j = x_f,j + C_f,xx * M^T (M C_f,xx M^T + C_ee)^-1 * (d_j-M*x_f,j)
% x_a,j = x_f,j + K * (d_j-M*x_f,j)

for iE = 1:EnKF.nE

    if sum([EnKF.Vel.loc, EnKF.Dir.loc, EnKF.TI.loc])>0
        % Calculate the distance of the OPs to the turbines
        OPs_tmp = ...
            EnKF.States_OP(:,EnKF.nStatesOP*(iE-1)+1:EnKF.nStatesOP*(iE-1)+3);
        distOPs = EnKF_distOPs(OPs_tmp,T.StartI);
    end
    
    % ========= Correct Velocity =========
    if EnKF.Vel.Correct
        % Calculate locaization
        LocCov = GaspariAndCohn1999(EnKF.Vel.cutOffLength,distOPs);

        % Set measurements and pollute with C_ee_Vel
        % d_Vel_j = d(:,1) + ...
        %     (randn(1,T.nT)*EnKF.Vel.C_ee_Vel_Chol)';
        
        % Kalman gain based on the mismatch in power generated
        K = EnKF_CalcKalmanGain('cxy_cyy_cee_loc', ...
            'Cxy',C_xy_Vel,'Cyy',C_yy_Vel,'Cee',EnKF.Vel.C_ee_Vel,...
            'Loc',LocCov);
        
        d_P_j = d_P + ...
            (randn(1,T.nT)*EnKF.Output.C_ee_Pow_Chol)';

        % Correction of the state
        EnKF.States.Vel(:,iE) = EnKF.States.Vel(:,iE) + ...
            K * (d_P_j - EnKF.Output.Pow(:,iE));% + ...
            %(randn(1,T.nT)*EnKF.Output.C_ee_Pow_Chol)');

    end

    % ========= Correct Direction =========
    if EnKF.Dir.Correct
        % Set measurements and pollute with C_ee_Dir
        d_Dir_j = d(:,2) + ...
            (randn(1,T.nT)*EnKF.Dir.C_ee_Dir_Chol)';

        % create localization covariance and multiply with state
        % error covariance matrix
        LocCov = GaspariAndCohn1999(EnKF.Dir.cutOffLength,distOPs);
        K = EnKF_CalcKalmanGain('c_cxx_cee_loc', ...
            'Cxx',C_xx_Dir,'Cee',EnKF.Dir.C_ee_Dir,...
            'Loc',LocCov,'C',C_phi);

        EnKF.States.Dir(:,iE) = EnKF.States.Dir(:,iE) + ...
            K * (d_Dir_j - EnKF.States.Dir(T.StartI,iE));

    end

    % ========= Correct TI =========
    if EnKF.TI.Correct
        % Set measurements
        d_TI_j = d(:,3) + ...
            (randn(1,T.nT)*EnKF.TI.C_ee_TI_Chol)';

        if EnKF.TI.loc
            % create localization covariance and multiply with state
            % error covariance matrix
            LocCov = GaspariAndCohn1999(EnKF.TI.cutOffLength,distOPs);

            K = (LocCov.*C_xx_TI)/(LocCov(T.StartI,:).*...
                C_xx_TI(T.StartI,:) + EnKF.TI.C_ee_TI);

            % Plot Localization influence
            if and(iS == 1, iE==1)
                EnKF_Vis = plotK_Localization(EnKF_Vis, T, Vis, ...
                    EnKF.TI.cutOffLength, 2);
            end
        else
            % No localization
            K = C_xx_TI/(C_xx_TI(T.StartI,:) + EnKF.TI.C_ee_TI);
        end

        % Plotting of the combined correction
        EnKF_Vis = plotK_combined(EnKF_Vis, T.nT, EnKF, Vis, K,...
            (d_TI_j - EnKF.States.TI(T.StartI,iE)),iE,3);

        % Correction
        EnKF.States.TI(:,iE) = EnKF.States.TI(:,iE) + ...
            K * (d_TI_j - EnKF.States.TI(T.StartI,iE));
    end

end

% Generate Plot
% CLC_Plot_EnsemblePower_vs_Measured(d_P, P,...
%     EnKF_calcPower(C_u, EnKF, T.D, paramFLORIS, T.StartI), [0,15])
% exportgraphics(gcf,[Sim.PathToSim filesep 'Results' ...
%             filesep 'EnKF201deg' filesep 'P_' num2str(time.current) '.png'])
% close(gcf)

% Clean up
clear K d_Dir_j d_TI_j d_Vel_j d C_xx_TI C_xx_Dir LocCov C_xy_Vel C_yy_Vel C_u C_phi d_P M