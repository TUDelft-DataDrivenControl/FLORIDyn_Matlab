function [J, M] = cost_function(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS)
%COST_FUNCTION Runs the FLORIDyn simulation and calculates a cost based on
%the simulation in- and output.

%% Run FLORIDyn simulation
[~,M,~] = FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);


%% Calculate the cost (decreasing is better)
J = -sum(M.("Power generated [MW]"));
end

