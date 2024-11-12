function T = correctDir(T,Wind,SimTime)
%CORRECTDIR Correction of the wind direction based on the OPs influencing
%the turbine.
%% Get Data
phi = getDataDir(Wind,T,SimTime);
%% Correct
for iT = 1:T.nT
    if isempty(T.dep{iT})
        T.States_WF(T.StartI(iT),2) = phi(iT);

        % OP Orientation = turbine wind direction
        if size(T.States_WF,2) == 4
            T.States_WF(T.StartI,4) = T.States_WF(T.StartI(iT),2);
        end
        continue
    end
    
    % Assign free wind speed based on what the influencing OPs carry
    if size(T.intOPs{iT},1)==1
        % Only one turbine influencing
        T.States_WF(T.StartI(iT),2) = ...
            T.States_WF(T.intOPs{iT}(1),2) * T.intOPs{iT}(2) + ...
            T.States_WF(T.intOPs{iT}(3),2) * T.intOPs{iT}(4);

    else
        % Weighted average of the wind speed state in the OPs that
        % influence this turbine. The weight is based on how much FLORIS
        % thinks this turbine is influenced by another wake. With the Gauss
        % model its currently based on the gaussian weight.
        if sum(T.weight{iT}) == 0
            % Turbine is actually not influenced -> skip
            T.States_WF(T.StartI(iT),2) = phi(iT);
        else
            % Sum the weighted wind speeds and the weights
            sum_wP = 0;
            sum_w  = 0;
            for iiT = 1:length(T.dep{iT})
                sum_wP = sum_wP + T.weight{iT}(iiT) * (...
                    T.States_WF(T.intOPs{iT}(iiT,1),2) .* T.intOPs{iT}(iiT,2) + ...
                    T.States_WF(T.intOPs{iT}(iiT,3),2) .* T.intOPs{iT}(iiT,4));
                
                sum_w = sum_w + T.weight{iT}(iiT);
            end
            % Divide to get weigthed average
            T.States_WF(T.StartI(iT),2) = sum_wP / sum_w;
        end
    end

    % OP Orientation = turbine wind direction
    if size(T.States_WF,2) == 4
        T.States_WF(T.StartI,4) = T.States_WF(T.StartI(iT),2);
    end
end


end

