function T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS)
% Initialize the simulation or load an initialized state
switch lower(Sim.Init)
    case 'init'
        
        if Sim.SaveInitState
            save([Sim.PathToSim 'T_init.mat'],'T')
        end
    case 'load'
        % load 
        try
            load([Sim.PathToSim 'T_init.mat'],'-mat','T')
        catch
            warning(['Could not load initT.mat from ' Sim.PathToSim ...
                '\nWill proceed with initialized data.'])
        end
end
end