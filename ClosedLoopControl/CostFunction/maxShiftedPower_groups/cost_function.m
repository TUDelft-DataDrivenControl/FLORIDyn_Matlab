function [J, M] = cost_function(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS)
%COST_FUNCTION Runs the FLORIDyn simulation and calculates a cost based on
%the simulation in- and output.

%% Special version!
% This version of the cost function is specifically designed for the 
% "yaw_time_shifted_group" controller

%% Run FLORIDyn simulation
[~,M,~] = FLORIDynCL(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

%% Calculate the cost (decreasing is better)
% Needs to be resorted to be [t1 k=0:ph,t2 k=0:ph, ...]'
p = reshape(M.("Power generated [MW]"),T.nT,[])';

% Shorten the relevant_t vector if not the entire ph was simulated
tmp = reshape(Con.shift.relevant_t,[],T.nT);
tmp = tmp(1:Sim.nSimSteps,:);

J = -tmp(:)' * p(:);
end

