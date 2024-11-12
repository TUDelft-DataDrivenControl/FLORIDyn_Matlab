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

