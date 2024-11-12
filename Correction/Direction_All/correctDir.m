function T = correctDir(T,Wind,SimTime)
%CORRECTDIR Correction of the wind direction
%% Get Data
phi = getDataDir(Wind,T,SimTime);
%% Correct
T.States_WF(:,2) = phi(1);

% OP Orientation = turbine wind direction
if size(T.States_WF,2) == 4
    T.States_WF(T.StartI,4) = phi(1);
end
end