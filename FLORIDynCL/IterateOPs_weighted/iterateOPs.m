function T = iterateOPs(T,Sim,paramFLORIS,paramFLORIDyn)
%ITERATEOPS propagates the OPs downstream and deletes the oldest ones.
% Weighted implementation
%   The idea is that OPs do propagate based on themselves but also on
%   other, neighbouring OPs. As a result, the free wind speed used for
%   propagation is a result of a weighted average function.

%% Save turbine OPs
tmpOPStates = T.States_OP(T.StartI,:);
tmpTStates  = T.States_T(T.StartI,:);
tmpWFSTates = T.States_WF(T.StartI,:);
%% Shift states
% Calculate downwind step and apply to wake coordinate system
WVel = getWeights(T, Sim.Dyn.Vel, Sim.TimeStep);
WDir = getWeights(T, Sim.Dyn.Dir, Sim.TimeStep);

WVel = WVel./sum(WVel,2);
WDir = WDir./sum(WDir,2);

% Determine advection speed based on advection model
propagationSpeed = WVel * T.States_WF(:,1);
turbineFreeSpeed = propagationSpeed(T.StartI);
switch Sim.Dyn.AdvectionModel
    case 'zong2020'
        % Advection speed based on U_adv(x) = 0.5*(U_inf + U_cen(x))
        propagationSpeed = propagationSpeed .*...
            getUadv(T.States_OP,T.States_T,T.States_WF,paramFLORIS,T.D(1));
    otherwise
        % Linear advection speed
        propagationSpeed = propagationSpeed * Sim.Dyn.Advection;
end


propagationDir   = WDir * T.States_WF(:,2); %<=== WIND DIR ADAPTATION


step_dw = Sim.TimeStep * propagationSpeed;
T.States_OP(:,4) = T.States_OP(:,4) + step_dw;

T.C_Vel = WVel(T.StartI,:);
T.C_Dir = WDir(T.StartI,:);

% Calculate crosswind step and apply to wake coordinate system
%   TODO If turbines with different diameters are used, this has to be run
%   individually (in parallel) for all turbines.
deflection  = Centerline(...
    T.States_OP,T.States_T,T.States_WF,paramFLORIS,T.D(1));
step_cw     = (deflection - T.States_OP(:,5:6))*1;

T.States_OP(:,5:6) = deflection*1;

% Apply dw and cw step to the world coordinate system
%phiW = angSOWFA2world(T.States_WF(:,2));
phiW = angSOWFA2world(propagationDir); %<=== WIND DIR ADAPTATION
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
T.States_WF(T.StartI,2) = propagationDir(T.StartI);
T.States_WF(T.StartI,1) = turbineFreeSpeed;
if size(T.States_WF,2) == 4
    % OP Orientation = turbine wind direction
    T.States_WF(T.StartI,4) = propagationDir(T.StartI);
end
%% Check if OPs are in order
for iT=1:T.nT
    [~,indOP] = sort(T.States_OP(T.StartI(iT)+(0:(T.nOP-1)),4));
    if ~issorted(indOP)
        orgInd = 0:(T.nOP-1);
        newInd = indOP - 1;
        
        % Find out where the new OPs are faster than the old ones
        newFaster =  T.States_WF(T.StartI(iT) + newInd,1)  > ...
            T.States_WF(T.StartI(iT)+orgInd,1);

%         figure(10)
%         hold on
%         scatter(T.States_OP(T.StartI(iT) + orgInd(newFaster),1),...
%             T.States_OP(T.StartI(iT) + orgInd(newFaster),2),20,...
%             T.States_WF(T.StartI(iT) + orgInd(newFaster),1),'filled')
%         hold off
%         axis equal
%         grid on
        
        %warning('still happening')
        T.States_OP(T.StartI(iT)+(0:(T.nOP-1)),:) = ...
            T.States_OP(T.StartI(iT) + indOP - 1,:);
        T.States_T(T.StartI(iT)+(0:(T.nOP-1)),:) = ...
            T.States_T(T.StartI(iT) + indOP - 1,:);
        T.States_WF(T.StartI(iT)+(0:(T.nOP-1)),:) = ...
            T.States_WF(T.StartI(iT) + indOP - 1,:);
    end
end

end

function W = getWeights(T, Sigma, timeStep)
% GETWEIGHTS generates a matrix with weights based on the distance of the
% OPs to each other. Wind direction is taken into account, also a temporal
% decay factor is implemented. Implementation inspired by Lejeune.
phiW = angSOWFA2world(T.States_WF(:,2));

distDW = ...
    cos(phiW).*(T.States_OP(:,1)-T.States_OP(:,1)') + ...
    sin(phiW).*(T.States_OP(:,2)-T.States_OP(:,2)');

distCW = ...
    sin(phiW).*(T.States_OP(:,1)-T.States_OP(:,1)') - ...
    cos(phiW).*(T.States_OP(:,2)-T.States_OP(:,2)');

% reduce calculation
% consider = or(abs(distDW)<3*Sigma.IterSigma_DW,...
%     abs(distCW)<3*Sigma.IterSigma_CW);
% 
% W = zeros(size(distDW));
% w1 = - distDW(consider).^2 / (2 * Sigma.IterSigma_DW^2);
% w2 = distCW(consider).^2 / (2 * Sigma.IterSigma_CW^2);
% W(consider) = exp(w1 - w2);

w1 = - distDW.^2 / (2 * Sigma.IterSigma_DW^2);
w2 = distCW.^2 / (2 * Sigma.IterSigma_CW^2);

W = exp(w1 - w2);

W = W .* repmat(...
    exp(-((0:T.nOP-1)*timeStep).^2 / (2 * Sigma.IterSigma_time^2)),...
    1,T.nT);
end

function W = getWeightsPar(T, Sigma, timeStep)
% GETWEIGHTS generates a matrix with weights based on the distance of the
% OPs to each other. Wind direction is taken into account, also a temporal
% decay factor is implemented. Implementation inspired by Lejeune.
phiW = angSOWFA2world(T.States_WF(:,2));

distDW = ...
    cos(phiW).*(T.States_OP(:,1)-T.States_OP(:,1)') + ...
    sin(phiW).*(T.States_OP(:,2)-T.States_OP(:,2)');

distCW = ...
    sin(phiW).*(T.States_OP(:,1)-T.States_OP(:,1)') - ...
    cos(phiW).*(T.States_OP(:,2)-T.States_OP(:,2)');

W = zeros(size(distCW(:)));
sig_dw = (2 * Sigma.IterSigma_DW^2);
sig_cw = (2 * Sigma.IterSigma_CW^2);
parfor i = 1:length(W)
    w1 = - distDW(i).^2 / sig_dw;
    w2 = distCW(i).^2 / sig_cw;
    
    W(i) = exp(w1 - w2);
end

W = reshape(W,size(distCW)) .* repmat(...
    exp(-((0:T.nOP-1)*timeStep).^2 / (2 * Sigma.IterSigma_time^2)),...
    1,T.nT);
end