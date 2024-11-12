function [Vis] = plotWeightedWindSpeed(T,Sim,Vis)
%PLOTAVERAGEDFLOWFIELD Summary of this function goes here
%   Detailed explanation goes here


fieldLims = Vis.FlowField.Lims;
fieldRes  = Vis.FlowField.Res;
xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);

[X,Y] = meshgrid(xAx,yAx);


%% Generate data
W = getWeights(X(:), Y(:), T, Sim);


%% Plotting
figure

mesh(X,Y,reshape(sum(W,2),size(X)),'FaceColor','interp','EdgeColor','none')
view(2)
axis equal
axis tight
colormap(inferno(16))


figure
W = W./sum(W,2);
mesh(X,Y,reshape(W * T.States_WF(:,1) * Sim.Dyn.Advection,size(X)),...
    'FaceColor','interp','EdgeColor','none')
view(2)
axis equal
axis tight
colormap(viridis(16))
colorbar
end

function W = getWeights(X, Y, T, Sim)
% GETWEIGHTS generates a matrix with weights based on the distance of the
% OPs to each other. Wind direction is taken into account, also a temporal
% decay factor is implemented. Implementation inspired by Lejeune.

distX = (X-T.States_OP(:,1)');
distY = (Y-T.States_OP(:,2)');
wPhi  = exp(- ...
    (distX.^2 + distY.^2) / (2 * Sim.Dyn.IterSigma_DW^2));
phi = (wPhi * T.States_WF(:,2))./sum(wPhi,2);


phiW = angSOWFA2world(phi);

distDW = ...
    cos(phiW).*(X-T.States_OP(:,1)') + ...
    sin(phiW).*(Y-T.States_OP(:,2)');

distCW = ...
    sin(phiW).*(X-T.States_OP(:,1)') - ...
    cos(phiW).*(Y-T.States_OP(:,2)');

W = exp(...
    - distDW.^2 / (2 * Sim.Dyn.IterSigma_DW^2) - ...
    distCW.^2 / (2 * Sim.Dyn.IterSigma_CW^2) );

W = W .* repmat(...
    exp(-((0:T.nOP-1)*Sim.TimeStep).^2 / (2 * Sim.Dyn.IterSigma_time^2)),...
    1,T.nT);
end