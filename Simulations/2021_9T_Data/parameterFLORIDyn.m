function paramFLORIDyn = parameterFLORIDyn()
%PARAMETERFLORIDYN Holds all parameters relevant to the FLORIDyn model
% Number of OPs per turbine
paramFLORIDyn.n_op = 200;

% Upstream area in which OPs are considered for interaction,
% normalized by the diameter of the turbine
%   Upstream distance (positive)
paramFLORIDyn.deltaUW = 10;
%   Downstream distance (positive)
paramFLORIDyn.deltaDW = .5;
%   Crosswind distance (3D -> 6D wide search area) (positive)
paramFLORIDyn.deltaCW = 3;

% Method for dynamic state change propagation
%   NOT IMPLEMENTED YET
paramFLORIDyn.dynStateChange = 'None';

% Temporary wind farm model
%   Two options exist to derive the TWF models; one is to set the OP 
%   orientation equal with the current wind direction. The other one is to
%   store the orientation as an additional state and derive the TWF based
%   on the initial particle orientation. The latter is in line with the
%   frozen turbulence theory, the former allows the full use of steady 
%   state models such as FLORIS to calculate the wake effect.
%   Settigs:
%   - 'homogeneous'     -> OP orientation = current wind direction
%   - 'heterogeneous'   -> OP orientation = wind direction at the time of
%                           the OP creation 
paramFLORIDyn.twf_model = 'heterogeneous';
end

