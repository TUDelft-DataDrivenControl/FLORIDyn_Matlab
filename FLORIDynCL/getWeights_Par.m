function W = getWeights_Par(X, Y, T, Dyn, TimeStep)
% GETWEIGHTS generates a matrix with weights based on the distance of the
% OPs to each other. Wind direction is taken into account, also a temporal
% decay factor is implemented. Implementation inspired by Lejeune.

distX = (X-T.States_OP(:,1)');
distY = (Y-T.States_OP(:,2)');
wPhi  = exp(- ...
    (distX.^2 + distY.^2) / (2 * Dyn.IterSigma_DW^2));
phi = (wPhi * T.States_WF(:,2))./sum(wPhi,2);


phiW = angSOWFA2world(phi);

distDW = ...
    cos(phiW).*(X-T.States_OP(:,1)') + ...
    sin(phiW).*(Y-T.States_OP(:,2)');

distCW = ...
    sin(phiW).*(X-T.States_OP(:,1)') - ...
    cos(phiW).*(Y-T.States_OP(:,2)');
Dyn_IterSigma_DW = 2*Dyn.IterSigma_DW^2;
Dyn_IterSigma_CW = 2*Dyn.IterSigma_CW^2;
W = zeros(size(distCW(:)));
parfor i = 1:length(W)
    W(i) = exp(...
        - distDW(i).^2 / Dyn_IterSigma_DW - ...
        distCW(i).^2 / Dyn_IterSigma_CW );
end

W = reshape(W,length(distX),[]) .* repmat(...
    exp(-((0:T.nOP-1)*TimeStep).^2 / (2 * Dyn.IterSigma_time^2)),...
    1,T.nT);
end
