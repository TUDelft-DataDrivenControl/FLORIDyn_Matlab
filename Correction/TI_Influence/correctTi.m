function T = correctTi(T,Wind,SimTime)
%CORRECTTI Correction of the turbulent intensity
%% Get Data
TI = getDataTI(Wind,T,SimTime);
%% Correct
for iT = 1:T.nT
    if isempty(T.dep{iT})
        T.States_WF(T.StartI(iT),3) = TI(iT);
        continue
    end
    
    % Assign amb. TI based on what the influencing OPs carry
    if size(T.intOPs{iT},1)==1
        % Only one turbine influencing
        T.States_WF(T.StartI(iT),3) = ...
            T.States_WF(T.intOPs{iT}(1),3) * T.intOPs{iT}(2) + ...
            T.States_WF(T.intOPs{iT}(3),3) * T.intOPs{iT}(4);
    else
        % More than one, currently only taking the mean, but could be more
        % sophisticated (e.g. only from uninfluenced wind turbines,
        % weighted by distance to OP)
        T.States_WF(T.StartI(iT),1) = mean(...
            T.States_WF(T.intOPs{iT}(:,1),3) .* T.intOPs{iT}(:,2) + ...
            T.States_WF(T.intOPs{iT}(:,3),3) .* T.intOPs{iT}(:,4));
    end
end
%T.States_WF(T.StartI,3) = TI;
end

