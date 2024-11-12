function [EnKF,trueOPpos] = EnKF_projectOntoTrueState(EnKF,Sim,T,stateIndex)
%EnKF_projectOntoTrueState projects the values of the 
switch stateIndex
    case 1
        % Velocity
        Dyn = Sim.Dyn.Vel;
        Val = EnKF.States.Vel;
    case 2
        % Wind direction
        Dyn = Sim.Dyn.Dir;
        Val = EnKF.States.Dir;
    case 3
        % Ambient turbulence intensity
        error('Amb. turb. not implemented yet!')
    otherwise
        error('Wrong state index, choose between 1,2 and 3.')
end
%% Get states and conversion matrices
valTrue = zeros(size(Val));
weightMat = cell(EnKF.nE,1);
%% Step 1: Calculate location of 'true' OPs
% trueOPpos = [...
%     mean(EnKF.States_OP(:,1:length(T.Names_OP):end),2),...
%     mean(EnKF.States_OP(:,2:length(T.Names_OP):end),2)];
%     
trueOPpos = ...
    reshape( ...
        mean(...
            reshape(...
                EnKF.States_OP,...
            length(T.Names_OP)*T.nT*T.nOP,[])...
        ,2),...
    T.nT*T.nOP,[]);
%% Step 2: Cycle through all ensembles and ...
% dist2tru = zeros(T.nOP*T.nT*EnKF.nE,1);
for iE = 1:EnKF.nE
    T = EnKF_AssignEnStates(EnKF,T,iE);

%% Step 2.1: Calculate the weights based on distance and OP age
    W = getWeightsT(trueOPpos(:,1), trueOPpos(:,2), T, Dyn, Sim.TimeStep);
    W = W./sum(W,2);
    weightMat{iE} = W;
%% Step 2.2: Calculate the values at the true OPs based on the weights
    valTrue(:,iE) = W * Val(:,iE);

    %% Test distance
%     distOPtoTrue = sqrt((trueOPpos(:,1)-T.States_OP(:,1)).^2 + ...
%         (trueOPpos(:,2)-T.States_OP(:,2)).^2);
%     dist2tru((iE-1)*T.nOP*T.nT+1:iE*T.nOP*T.nT) = distOPtoTrue;
end



%% Store in EnKF struct
switch stateIndex
    case 1
        % Velocity
        EnKF.Vel.W = weightMat;
        EnKF.States.Vel = valTrue;

        %% PLOTTING OF DISTANCES
        
%         figure(10)
%         hold on

%         binWidth = .1;
%         [N,edges] = histcounts(dist2tru./32.8,0:binWidth:2,'Normalization','cdf');
%         if(~isfield(EnKF,'NormDistAll'))
%             EnKF.NormDistAll = N;
%         else
%             EnKF.NormDistAll = [N; EnKF.NormDistAll];
%         end

%         plot(edges(1:end-1)+binWidth/2, N, '--*', 'LineWidth',1)
%         histogram(dist2tru./32.8000,25,'Normalization','probability',...
%             'DisplayStyle','stairs','LineWidth',1)
%         grid on
%         xlim([0,2])
%         ylim([0,1])
%         xlabel('Normalized distance')
%         ylabel('Cumalative density function estimate')
%         hold off
%         title('Entire simulation space')
% 
%         figure(11)
%         hold on

%         withinSim = and(trueOPpos(:,1)<3000,(trueOPpos(:,2)<3000));
%         [N,edges] = histcounts(dist2tru(withinSim)./32.8,0:binWidth:2,...
%             'Normalization','cdf');
%         if(~isfield(EnKF,'NormDistRed'))
%             EnKF.NormDistRed = N;
%         else
%             EnKF.NormDistRed = [N; EnKF.NormDistRed];
%         end

%         plot(edges(1:end-1)+binWidth/2, N, '--*', 'LineWidth',1)
%         grid on
%         xlim([0,2])
%         ylim([0,1])
%         xlabel('Normalized distance')
%         ylabel('Cumalative density function estimate')
%         hold off
%         title('Relevant simulation space')
    case 2
        % Wind direction
        EnKF.Dir.W = weightMat;
        EnKF.States.Dir = valTrue;
    case 3
        % Ambient turbulence intensity
        error('Amb. turb. not implemented yet!')
    otherwise
        error('Wrong state index, choose between 1,2 and 3.')
end
end

%% Additional functions
function W = getWeights(X, Y, T, Dyn, TimeStep)
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

W = exp(...
    - distDW.^2 / (2 * Dyn.IterSigma_DW^2) - ...
    distCW.^2 / (2 * Dyn.IterSigma_CW^2) );

W = W .* repmat(...
    exp(-((0:T.nOP-1)*TimeStep).^2 / (2 * Dyn.IterSigma_time^2)),...
    1,T.nT);
end