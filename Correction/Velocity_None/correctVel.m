function [T,Wind] = correctVel(T,Wind,SimTime,paramFLORIS,tmpM)
%CORRECTVEL Correction of the Velocity state of the OPs

%% get data
[U, Wind] = getDataVel(Wind,T,SimTime,tmpM,paramFLORIS);
%% Correct
T.States_WF(T.StartI,1) = U;
end