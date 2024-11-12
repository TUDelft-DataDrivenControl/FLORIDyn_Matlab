%% Control script 
% Create true state from EnKF data

% Optimize
[Con, CLC] = controller(T, Wind, Sim, Con, Vis, paramFLORIDyn,...
    paramFLORIS,CLC, Sim.StartTime);