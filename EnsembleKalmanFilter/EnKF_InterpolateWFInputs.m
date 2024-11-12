function [EnKF,Wind] = EnKF_InterpolateWFInputs(EnKF,Wind,Sim,T)
%ENKF_INTERPOLATEWFINPUTS Method to write Wind.Vel .Dir .TI in such a way
%that the values get interpolated

% Save current state for state based input methods
%   NOTE: This only works if all ensembles have the same state, otherwise
%   this just uses the state from the last ensemble
U   = T.States_WF(T.StartI,1);
phi = T.States_WF(T.StartI,2);

for iW = 1:EnKF.nS+1
    SimTime = Sim.StartTime+(iW-1)*Sim.TimeStep;
    Wind.Vel(iW,1) = SimTime;
    Wind.Dir(iW,1) = SimTime;
    Wind.TI(iW,1) = SimTime;
    %% Wind speed
    
    if strcmp(EnKF.Wind.Input_Vel,'I_and_I')
        [U,EnKF.Wind.Vel] = getWindSpeedT(EnKF.Wind.Vel,(1:T.nT)',SimTime,...
            T.States_WF(T.StartI,2),paramFLORIS.p_p);
        if (SimTime-EnKF.Wind.Vel.StartTime)>EnKF.Wind.Vel.WSE.Offset
            % Ufree = Ueff/reduction
            %U = U./tmpM(:,1);
            warning(['The EnKF can not yet use the I&I as intended, '...
                'proceeds with U_eff=U_free'])
        end
    elseif strcmp(EnKF.Wind.Input_Vel,'ZOH_wErrorCov')
        U = getWindSpeedT(U, EnKF.Wind.Vel.ColSig);
    elseif strcmp(EnKF.Wind.Input_Vel,'RW_with_Mean')
        U = getWindSpeedT(U,EnKF.Wind.Vel);
    else
        U = getWindSpeedT(EnKF.Wind.Vel,(1:T.nT)',SimTime);
    end
    Wind.Vel(iW,2:end) = U';
    
%     Wind.Vel(iW,2:end) = getWindSpeedT(EnKF.Wind.Vel,1:T.nT,...
%         SimTime)';
    
    %% Wind direction
    if strcmp(EnKF.Wind.Input_Vel,'RW_with_Mean')
        phi = getWindDirT(phi,EnKF.Wind.Dir);
    else
        phi = getWindDirT(EnKF.Wind.Dir,(1:T.nT)',SimTime);
    end
    
    Wind.Dir(iW,2:end) = phi';
%     Wind.Dir(iW,2:end) = getWindDirT(EnKF.Wind.Dir,1:T.nT,...
%         SimTime)';
    
    %% Ambient turbulence intensity
    Wind.TI(iW,2:end) = getWindTiT(EnKF.Wind.TI,1:T.nT,...
        SimTime)';
end
end

