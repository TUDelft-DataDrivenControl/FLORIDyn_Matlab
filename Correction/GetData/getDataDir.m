function phi = getDataDir(Wind,T,SimTime)
%GETDATADIR Reads wind data and returns the current phi for all turbines

switch Wind.Input_Dir
    case 'RW_with_Mean'
        phi = getWindDirT(T.States_WF(T.StartI,2),Wind.Dir);
    case 'EnKF_ZOH'
        phi = T.States_WF(T.StartI,2);
    case 'EnKF_RW'
        phi = T.States_WF(T.StartI,2);
        phi = phi + 0*(randn(1,length(phi))*Wind.Dir.CholSig)';
    case 'EnKF_InterpTurbine'
        phi = getWindDirT_EnKF(Wind.Dir,(1:T.nT)',SimTime);
    case 'CLC_weighted_ZOH'
        % Apply weighted states on the turbine
        phi = T.C_Dir*T.States_WF(:,2);
    otherwise
        phi = getWindDirT(Wind.Dir,(1:T.nT)',SimTime);
end

end

