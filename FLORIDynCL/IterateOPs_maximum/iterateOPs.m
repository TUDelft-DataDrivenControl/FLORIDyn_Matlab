function T = iterateOPs(T,Sim,paramFLORIS,paramFLORIDyn)
%ITERATEOPS propagates the OPs downstream and deletes the oldest ones.


%% Save turbine OPs
tmpOPStates = T.States_OP(T.StartI,:);
tmpTStates  = T.States_T(T.StartI,:);
tmpWFSTates = T.States_WF(T.StartI,:);
%% Shift states
% Calculate downwind step and apply to wake coordinate system
step_dw = Sim.TimeStep * T.States_WF(:,1) * Sim.Dyn.Advection;
T.States_OP(:,4) = T.States_OP(:,4) + step_dw;

% Calculate crosswind step and apply to wake coordinate system
%   TODO If turbines with different diameters are used, this has to be run
%   individually (in parallel) for all turbines.
deflection  = Centerline(...
    T.States_OP,T.States_T,T.States_WF,paramFLORIS,T.D(1));
step_cw     = (deflection - T.States_OP(:,5:6))*1;

T.States_OP(:,5:6) = deflection*1;

% Apply dw and cw step to the world coordinate system
phiW = angSOWFA2world(T.States_WF(:,2));

T.States_OP(:,1) = T.States_OP(:,1) + ...
    cos(phiW) .* step_dw - ...
    sin(phiW) .* step_cw(:,1);
T.States_OP(:,2) = T.States_OP(:,2) + ...
    sin(phiW) .* step_dw + ...
    cos(phiW) .* step_cw(:,1);
T.States_OP(:,3) = T.States_OP(:,3) + step_cw(:,2);

%% Circshift & Init first OPs
%   OPs
T.States_OP = circshift(T.States_OP,1);
T.States_OP(T.StartI,:) = tmpOPStates;
%   Turbines
T.States_T  = circshift(T.States_T,1);
T.States_T(T.StartI,:)  = tmpTStates;
%   Wind Farm
T.States_WF = circshift(T.States_WF,1);
T.States_WF(T.StartI,:) = tmpWFSTates;

%% Check if OPs are in order and overwrite the slower one
for iT=1:T.nT
    [~,indOP] = sort(T.States_OP(T.StartI(iT)+(0:(T.nOP-1)),4));
    if ~issorted(indOP)
        orgInd = 0:(T.nOP-1);
        newInd = indOP - 1;
        
        % Find out where the new OPs are faster than the old ones
        newFaster =  T.States_WF(T.StartI(iT) + newInd,1)  > ...
            T.States_WF(T.StartI(iT)+orgInd,1);
        
        % Overwrite states of slower OP
        T.States_OP(T.StartI(iT) + orgInd(newFaster), :) = ...
            T.States_OP(T.StartI(iT) + newInd(newFaster),:);
        T.States_T(T.StartI(iT) + orgInd(newFaster), :) = ...
            T.States_T(T.StartI(iT) + newInd(newFaster),:);
        T.States_WF(T.StartI(iT) + orgInd(newFaster), :) = ...
            T.States_WF(T.StartI(iT) + newInd(newFaster),:);
    end
end
end

