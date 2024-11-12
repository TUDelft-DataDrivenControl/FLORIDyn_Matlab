function EnKF = EnKF_projectOntoEnsembleState(EnKF,stateIndex)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for iE = 1:EnKF.nE
    switch stateIndex
        case 1
            %EnKF.States.Vel(:,iE) = EnKF.Vel.W{iE}\EnKF.States.Vel(:,iE);
            EnKF.States.Vel(:,iE) = ...
                linsolve(EnKF.Vel.W{iE},EnKF.States.Vel(:,iE));
        case 2
            EnKF.States.Dir(:,iE) = EnKF.Dir.W{iE}\EnKF.States.Dir(:,iE);
        case 3
            error('Not implemented yet')
    end
end
end