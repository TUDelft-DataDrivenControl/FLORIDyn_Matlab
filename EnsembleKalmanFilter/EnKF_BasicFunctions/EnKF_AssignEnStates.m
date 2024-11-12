function T = EnKF_AssignEnStates(EnKF,T,iE)
%ENKF_ASSIGNENSTATES Assign States if the wind field state is being
%corrected for the next simulation run.

%% Ensemble States
if EnKF.Vel.Correct
    T.States_WF(:,1) = EnKF.States.Vel(:,iE);
end
if EnKF.Dir.Correct
    T.States_WF(:,2) = EnKF.States.Dir(:,iE);
end
if EnKF.TI.Correct
    T.States_WF(:,3) = EnKF.States.TI(:,iE);
end

%% OP States
T.States_OP = ...
    EnKF.States_OP(:,EnKF.nStatesOP*(iE-1)+1:EnKF.nStatesOP*iE);

%% T States
T.States_T = ...
    EnKF.States_T(:,EnKF.nStatesT*(iE-1)+1:EnKF.nStatesT*iE);
end

% === TODO ===
% Uncorrected states are just passed on now, from one ensemble to the next.
% This mostly concerns the uncorrected WF states (e.g. Ambient TI)
%   No issue if the value does not change anyways.