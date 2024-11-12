% CLC EnKF script
tic
% Set inputs
Sim.StartTime   = time.current - time.stepEnKF;
Sim.EndTime     = time.current;
Sim.nSimSteps   = EnKF.nS;

% Prepare parallel computing
EnKFStatesVel = EnKF.States.Vel;
EnKFStatesDir = EnKF.States.Dir;
EnKFStates_OP = reshape(EnKF.States_OP,[],EnKF.nE);
EnKFStates_T  = reshape(EnKF.States_T,[],EnKF.nE);
EnKFnStatesOP = EnKF.nStatesOP;
EnKFnStatesT  = EnKF.nStatesT;
EnKFOutputPow = EnKF.Output.Pow;
EnKFM         = zeros(Sim.nSimSteps * T.nT * 6,EnKF.nE);

parfor iE = 1:EnKF.nE
    % ===== Assign relevant ensemble states
    T_tmp = T;

    % Ensemble States
    T_tmp.States_WF(:,1) = EnKFStatesVel(:,iE);
    T_tmp.States_WF(:,2) = EnKFStatesDir(:,iE);
    
    % OP States
    T_tmp.States_OP = reshape(EnKFStates_OP(:,iE),[],EnKFnStatesOP);
        %EnKF.States_OP(:,EnKF.nStatesOP*(iE-1)+1:EnKF.nStatesOP*iE);
    
    % T States
    T_tmp.States_T = reshape(EnKFStates_T(:,iE),[],EnKFnStatesT);
        %EnKF.States_T(:,EnKF.nStatesT*(iE-1)+1:EnKF.nStatesT*iE);

    %T_tmp = EnKF_AssignEnStates(EnKF,T,iE);
    
    % =========== Run simulation ================
    [T_tmp,M,~] = ...
        FLORIDynCL(T_tmp,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
    
    % ===== Stores relevant measurements and states
    % Store Ensemble states
    EnKFStatesVel(:,iE) = T_tmp.States_WF(:,1);
    EnKFStatesDir(:,iE) = T_tmp.States_WF(:,2);
    
    % Store OP states
    EnKFStates_OP(:,iE) = T_tmp.States_OP(:);
    
    % Store T states
    EnKFStates_T(:,iE)  = T_tmp.States_T(:);
    
    % Store measurements
    EnKFM(:,iE) = reshape(table2array(M),[],1);
    %EnKF.M{iE} = [EnKF.M{iE}; M];
    EnKFOutputPow(:,iE) = ...
        table2array(M(end-T.nT+1:end,6));
    %EnKF = EnKF_StoreEnStates(EnKF,M, T_tmp, iE);
end

% Limit
if EnKF.Vel.Limits
    EnKFStatesVel(EnKFStatesVel<EnKF.Vel.Min) = EnKF.Vel.Min;
    EnKFStatesVel(EnKFStatesVel>EnKF.Vel.Max) = EnKF.Vel.Max;
end

% Collect parallel computing
EnKF.States.Vel = EnKFStatesVel;
EnKF.States.Dir = EnKFStatesDir;
EnKF.States_OP  = reshape(EnKFStates_OP,[], EnKF.nE * EnKF.nStatesOP);
EnKF.States_T   = reshape(EnKFStates_T, [], EnKF.nE * EnKF.nStatesT);
EnKF.Output.Pow = EnKFOutputPow;
for iE = 1:EnKF.nE
    EnKF.M{iE} = [EnKF.M{iE}; ...
        array2table(reshape(EnKFM(:,iE),Sim.nSimSteps * T.nT,[]),...
        'VariableNames',...
        {'Time [s]','Foreign Reduction [%]','Added turbulence [%]',...
        'Effective Wind Speed [ms^-1]','Free Wind Speed [ms^-1]', ...
        'Power generated [MW]'})];
end
time_simulate = toc;
clear EnKFStatesVel EnKFStatesDir EnKFStates_OP EnKFStates_T EnKFM
%% Combination
tic
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
time_combine = toc;
%% Get measurements from validation source
%   Calls the functions to get wind speed, direction and amb.
%   turbulence as the normal simulation would and gets [U,phi,TI]
[d, EnKF] = EnKF_GetWFInputs(EnKF,Sim,T,paramFLORIS);
d_phi = d(:,2);
% Calculate reduction of free wind speed at the turbines & apply to meas.
% r = EnKF_GetAveForeignReduction(EnKF.M,T.nT); 
% d(:,1) = d(:,1)./(r*10^(-2)); 

if powSOWFA == 0
    % Running closed loop
    %   Retrieve power generated from ZeroMQ message
    %   (Air density error correction not needed!)
    d_P   = generatorPowerArray'./(10^6);

    %   Retrieve wind direction from probe data
    if startsWith(EnKF.Dir.MeasurementSource, 'CLC Probes ')
        try
            d_phi = CLC_read_probe(EnKF.Dir.PathToProbeFile, ...
                EnKF.Dir.MeasurementAverageTime, EnKF.Dir.MeasurementDefault)';
        catch
            disp('ERROR READING WIND DIRECTION PROBE FILE!')
            d_phi = ones(T.nT,1)*EnKF.Dir.MeasurementDefault;
        end
        if endsWith(EnKF.Dir.MeasurementSource, 'Averaged')
            d_phi(:) = mean(d_phi);
        end
    end
else
    % Running offline
    d_P = interp1(powSOWFA(1:T.nT:end,2),...
        reshape(powSOWFA(:,3),T.nT,[])',Sim.EndTime)';
end



%% Correct Ensemble states
% [1] Eq.4.37
% x_a,j = x_f,j + C_f,xx * M^T (M C_f,xx M^T + C_ee)^-1 * (d_j-M*x_f,j)
% x_a,j = x_f,j + K * (d_j-M*x_f,j)
tic
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

        % EnKF.Output.Pow IS NOT SET!!! IS ZERO!!!

    end

    % ========= Correct Direction =========
    if EnKF.Dir.Correct
        % Set measurements and pollute with C_ee_Dir
        d_Dir_j = d_phi + ...
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
time_correct = toc;

%% Store data
if EnKF.StoreData
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'EnKF_Power'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end

    % Store Power
    writematrix(d_P,[folder_name filesep...
        pad(num2str(time.current),7,"left",'0') '_power_measured.txt'])
    writematrix(P,[folder_name filesep...
        pad(num2str(time.current),7,"left",'0') '_power_redicted.txt'])
    writematrix(EnKF_calcPower(C_u, EnKF, T.D, paramFLORIS, T.StartI),...
        [folder_name filesep...
        pad(num2str(time.current),7,"left",'0') '_power_corrected.txt'])
    
   % Store Wind direction
    writematrix(d_phi,[folder_name filesep...
        pad(num2str(time.current),7,"left",'0') '_dir_measured.txt'])

    % Store time
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'EnKF_Time'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end
    writematrix([time.current time_simulate time_combine time_correct],...
        [folder_name filesep 'time_sim_com_cor.csv'],"WriteMode","append")
    
end

% Clean up
clear K d_Dir_j d_TI_j d_Vel_j d C_xx_TI C_xx_Dir LocCov C_xy_Vel C_yy_Vel C_u C_phi d_P M time_simulate time_combine time_correct d_phi