function T = correctTi(T,Wind,SimTime)
%CORRECTTI Correction of the turbulent intensity
%% Get Data
TI = getDataTI(Wind,T,SimTime);
%% Correct
T.States_WF(T.StartI,3) = TI;
end

