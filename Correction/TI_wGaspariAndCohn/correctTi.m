function T = correctTi(T,Wind,SimTime)
%CORRECTTI Correction of the turbulent intensity
warning(['correctTi() is not yet correctly implemented for Graspari'...
    ' and Cohn influence correction'])
%% Get Data
TI = getDataTI(Wind,T,SimTime);
%% Correct
T.States_WF(T.StartI,3) = TI;
end

