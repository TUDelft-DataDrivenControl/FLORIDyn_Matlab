function [T,Wind] = correctVel(T,Wind,SimTime,paramFLORIS,tmpM)
%CORRECTVEL Correction of the Velocity state of the OPs

%% get data
[U, Wind] = getDataVel(Wind,T,SimTime,tmpM,paramFLORIS);
%% Correct
for iT = 1:T.nT
    if isempty(T.dep{iT})
        T.States_WF(T.StartI(iT),1) = U(iT);
        continue
    end
    
    % Assign free wind speed based on what the influencing OPs carry
    if size(T.intOPs{iT},1)==1
        % Only one turbine influencing
        T.States_WF(T.StartI(iT),1) = ...
            T.States_WF(T.intOPs{iT}(1),1) * T.intOPs{iT}(2) + ...
            T.States_WF(T.intOPs{iT}(3),1) * T.intOPs{iT}(4);
    else
        % Weighted average of the wind speed state in the OPs that
        % influence this turbine. The weight is based on how much FLORIS
        % thinks this turbine is influenced by another wake. With the Gauss
        % model its currently based on the gaussian weight.
        if sum(T.weight{iT}) == 0
            % Turbine is actually not influenced -> skip
            T.States_WF(T.StartI(iT),1) = U(iT);
            continue
        else
            % Sum the weighted wind speeds and the weights
            sum_wU = 0;
            sum_w  = 0;
            for iiT = 1:length(T.dep{iT})
                sum_wU = sum_wU + T.weight{iT}(iiT) * (...
                    T.States_WF(T.intOPs{iT}(iiT,1),1) .* T.intOPs{iT}(iiT,2) + ...
                    T.States_WF(T.intOPs{iT}(iiT,3),1) .* T.intOPs{iT}(iiT,4));
                
                sum_w = sum_w + T.weight{iT}(iiT);
            end
            % Divide to get weigthed average
            T.States_WF(T.StartI(iT),1) = sum_wU / sum_w;
        end
    end
end
end

