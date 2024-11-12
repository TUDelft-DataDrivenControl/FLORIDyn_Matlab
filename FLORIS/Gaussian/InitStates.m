function [States_OP, States_T, States_WF] = InitStates(T,Wind,InitTurb,paramFLORIS,Sim)
%% Fill the states with initial values and locate OPs at a meaningful pos.
States_OP   = T.States_OP;
States_T    = T.States_T;
States_WF   = T.States_WF;
nT          = T.nT;
nOP         = T.nOP;
deltaT      = Sim.TimeStep;
startTime   = Sim.StartTime;

for iT = 1:nT
    % Retrieve wind field data
    if strcmp(Wind.Input_Vel,'I_and_I')
        [U,~] = getWindSpeedT(Wind.Vel,iT,startTime);
    elseif or(strcmp(Wind.Input_Vel,'ZOH_wErrorCov'),...
            strcmp(Wind.Input_Vel,'RW_with_Mean'))
        U = Wind.Vel.Init;
    else
        U     = getWindSpeedT(Wind.Vel,iT,startTime);
    end
    
    if strcmp(Wind.Input_Dir,'RW_with_Mean')
        phiS = Wind.Dir.Init;
    else
        phiS = getWindDirT(Wind.Dir,iT,startTime);
    end
    
    TI   = getWindTiT(Wind.TI,iT,startTime);
    
    rangeOPs = (iT-1)*nOP+1:iT*nOP;
    %Iinitialize the States of the OPs and turbines
    States_WF(rangeOPs,1) = U;
    States_WF(rangeOPs,2) = phiS; % SOWFA angle
    States_WF(rangeOPs,3) = TI;
    % Add orientation if used
    if length(T.Names_WF) == 4
        States_WF(rangeOPs,4) = phiS;
    end
    
    % Downwind distance (wake coord)
    States_OP(rangeOPs,4) = (0:nOP-1)' * deltaT * U;
    
    % Init turbine states
    States_T(rangeOPs,:) = ones(nOP,1)*InitTurb(iT,:);
    
    % Crosswind position
    States_OP(rangeOPs,5:6) = ...
        Centerline(States_OP(rangeOPs,:),States_T(rangeOPs,:),...
            States_WF(rangeOPs,:),paramFLORIS,T.D(iT));
    
    % Convert wind dir in fitting radians
    phiW = angSOWFA2world(States_WF(rangeOPs,2));
    
    % World coordinate position x0 and y0 including tower base and nac pos
    States_OP(rangeOPs,1) = ...
        cos(phiW) .* States_OP(rangeOPs,4) - ...
        sin(phiW) .* States_OP(rangeOPs,5) + ...
        T.posBase(iT,1) + T.posNac(iT,1);
    
    States_OP(rangeOPs,2) = ...
        sin(phiW) .* States_OP(rangeOPs,4) + ...
        cos(phiW) .* States_OP(rangeOPs,5) + ...
        T.posBase(iT,2) + T.posNac(iT,2);
    
    States_OP(rangeOPs,3) = States_OP(rangeOPs,6) + ...
                            T.posBase(iT,3) + T.posNac(iT,3);
end

end

