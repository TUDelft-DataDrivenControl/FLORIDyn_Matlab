function TI = getDataTI(Wind,T,SimTime)
%GETDATATI retrieves the data for the ambient turbulence intensity

switch Wind.Input_TI
    case 'EnKF_InterpTurbine'
        TI = getWindTiT_EnKF(Wind.TI,(1:T.nT)',SimTime);
    case 'EnKF_ZOH'
        TI = T.States_WF(T.StartI,3);
    case 'EnKF_RW'
        TI = T.States_WF(T.StartI,3);
        TI = TI + (randn(1,length(phi))*Wind.TI.CholSig)';
    case 'CLC_weighted_ZOH'
        % Apply weighted states on the turbine
        TI = T.C_TI*T.States_WF(:,3);
    otherwise
        TI = getWindTiT(Wind.TI,(1:T.nT)',SimTime);
end
end

