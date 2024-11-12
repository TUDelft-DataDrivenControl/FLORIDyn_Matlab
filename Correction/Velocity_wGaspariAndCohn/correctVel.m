function [T,Wind] = correctVel(T,Wind,SimTime,paramFLORIS,tmpM)
%CORRECTVEL Correction of the Velocity state of the OPs
warning(['correctVel() is not yet correctly implemented for Graspari'...
    ' and Cohn influence correction'])
%% get data
[U, Wind] = getDataVel(Wind,T,SimTime,tmpM,paramFLORIS);
%% Correct
T.States_WF(T.StartI,1) = U;
end

