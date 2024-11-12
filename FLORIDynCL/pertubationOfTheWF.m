function T = pertubationOfTheWF(T,Wind)
%% pertubationOfTheWF adds noise to the entire wind field state
% Velocity
if Wind.Pertubation.Vel
    T.States_WF(:,1) = T.States_WF(:,1) + ...
        Wind.Pertubation.VelSigma * randn(T.nOP*T.nT,1);
end
% Direction
if Wind.Pertubation.Dir
    T.States_WF(:,2) = T.States_WF(:,2) + ...
        Wind.Pertubation.DirSigma * randn(T.nOP*T.nT,1);
end
% Turbulence intensity
if Wind.Pertubation.TI
    T.States_WF(:,3) = T.States_WF(:,3) + ...
        Wind.Pertubation.TISigma * randn(T.nOP*T.nT,1);
end

end