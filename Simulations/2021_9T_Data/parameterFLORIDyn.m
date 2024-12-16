function paramFLORIDyn = parameterFLORIDyn()
%PARAMETERFLORIDYN Holds all parameters relevant to the FLORIDyn model
% Number of OPs per turbine
paramFLORIDyn.n_op = 130;

% Upstream area in which turbines are considered for interaction,
% normalized by the diameter of the turbine
%   Upstream distance (positive)
paramFLORIDyn.deltaDW = 10;
%   Crosswind distance (3D -> 6D wide search area) (positive)
paramFLORIDyn.deltaCW = 3;

% Method for dynamic state change propagation
%   NOT IMPLEMENTED YET
paramFLORIDyn.dynStateChange = 'None';
end

