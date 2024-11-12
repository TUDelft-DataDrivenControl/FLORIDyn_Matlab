function W = getWeights(X, Y, T, Dyn, TimeStep)
% GETWEIGHTS generates a matrix with weights based on the distance of the
% OPs to each other. Wind direction is taken into account, also a temporal
% decay factor is implemented. Implementation inspired by Lejeune.

max_sig = max(Dyn.IterSigma_DW,Dyn.IterSigma_CW);          % Addition to reduce calculations

distX = (X-T.States_OP(:,1)');
distY = (Y-T.States_OP(:,2)');

consider = or(abs(distX)<3*max_sig, abs(distY)<3*max_sig); % Addition to reduce calculations

wPhi = zeros(distX);
% wPhi  = exp(- ...
%     (distX.^2 + distY.^2) / (2 * Dyn.IterSigma_DW^2));

wPhi(consider)  = exp(- ...
    (distX(consider).^2 + distY(consider).^2) / (2 * Dyn.IterSigma_DW^2));

phi = (wPhi * T.States_WF(:,2))./sum(wPhi,2);
                                    
phiW = angSOWFA2world(phi);

% distDW = ...
%     cos(phiW).*(X-T.States_OP(:,1)') + ...
%     sin(phiW).*(Y-T.States_OP(:,2)');
% 
% distCW = ...
%     sin(phiW).*(X-T.States_OP(:,1)') - ...
%     cos(phiW).*(Y-T.States_OP(:,2)');



distDW = ...
    cos(phiW).*distX(consider) + ...
    sin(phiW).*distY(consider);

distCW = ...
    sin(phiW).*distX(consider) - ...
    cos(phiW).*distY(consider);

W = zeros(size(distX));

W(consider) = exp(...
    - distDW.^2 / (2 * Dyn.IterSigma_DW^2) - ...
    distCW.^2 / (2 * Dyn.IterSigma_CW^2) );

% W = exp(...
%     - distDW.^2 / (2 * Dyn.IterSigma_DW^2) - ...
%     distCW.^2 / (2 * Dyn.IterSigma_CW^2) );

W = W .* repmat(...
    exp(-((0:T.nOP-1)*TimeStep).^2 / (2 * Dyn.IterSigma_time^2)),...
    1,T.nT);
end
