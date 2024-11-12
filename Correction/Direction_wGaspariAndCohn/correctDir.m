function T = correctDir(T,Wind,SimTime)
%CORRECTDIR Correction of the wind direction
warning(['correctDir() is not yet correctly implemented for Graspari'...
    ' and Cohn influence correction'])
%% Get Data
phi = getDataDir(Wind,T,SimTime);
%% Correct
T.States_WF(T.StartI,2) = phi;

% OP Orientation = turbine wind direction
if size(T.States_WF,2) == 4
    T.States_WF(T.StartI,4) = T.States_WF(T.StartI(iT),2);
end
end

