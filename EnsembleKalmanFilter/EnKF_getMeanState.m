function [T, M] = EnKF_getMeanState(EnKF, T, paramFLORIDyn)
%ENKF_GETMEANSTATE Returns the current estimated state of the system

%% Ensemble States
if EnKF.Vel.Correct
    T.States_WF(:,1) = mean(EnKF.States.Vel,2);
    T.States_Var_Vel = var(EnKF.States.Vel,[],2);
end
if EnKF.Dir.Correct
    T.States_WF(:,2) = mean(EnKF.States.Dir,2);
    T.States_Var_Dir = var(EnKF.States.Dir,[],2);
end
if EnKF.TI.Correct
    T.States_WF(:,3) = mean(EnKF.States.TI,2);
end

%% OP States
for iS = 1:EnKF.nStatesOP
    T.States_OP(:,iS) = ...
        mean(EnKF.States_OP(:,iS:EnKF.nStatesOP:end),2);
end

%% T States
for iS = 1:EnKF.nStatesT
    T.States_T(:,iS) = ...
        mean(EnKF.States_T(:,iS:EnKF.nStatesOP:end),2);
end

%% Calculate dependencies between turbines
%   ======== Find groups
T.dep = findTurbineGroups(T, paramFLORIDyn);

% %   ======= Interpolate OPs
% T.intOPs = interpolateOPs(T);

%% NEEDS TO BE FIXED TO... MEAN?
M = EnKF.M{1};
end

