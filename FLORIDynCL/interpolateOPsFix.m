function intOPs = interpolateOPsFix(T)
% ONLY WORKS WITH A WIND DIRECTION OF 270 DEG!

%INTERPOALTEOPS Interpolates the OPs at the closest location of the
%turbines which are then used to set up the temporary wind farms.

intOPs = cell(T.nT,1);
StartI = T.StartI;

for iT = 1:T.nT % Every turbine
    % [id OP1, rat OP1, id OP2, rat OP2]
    intOPs{iT} = zeros(length(T.dep{iT}),4);
    
    for iiT = 1:length(T.dep{iT}) % every influencing turbine
        % actual turbine index
        iiaT = T.dep{iT}(iiT);
        
        % distances of influencing OPs to turbine
        dist = sqrt(sum(...
            (T.States_OP(StartI(iiaT)+(0:T.nOP-1),1)... % <= only dep. on X
            - T.posBase(iT,1)).^2,2));                  % <= only dep. on X
        
        % get closest OP
        [~,i_min] = sort(dist);
        
        % get second closest
        if i_min(1) == 1
            % first OP (unlikely)
            intOPs{iT}(iiT,:) = [StartI(iiaT), 1, StartI(iiaT)+1, 0];
        elseif i_min(1) == T.nOP
            % last OP (possible)
            intOPs{iT}(iiT,:) = [StartI(iiaT)+T.nOP-2, 0, ...
                StartI(iiaT)+T.nOP-1, 1];
        else
            %% use two closest OPs to interpolate influence
            indOP1 = StartI(iiaT) - 1 + i_min(1);
            indOP2 = StartI(iiaT) - 1 + i_min(2);
            % Relative Distance: d = -((b-a)⊤(a-c))/((b-a)⊤(b-a))
            % Interpolation location: a+(b-a)d
            % a: OP close, b: OP far, c: turb_location
            % Interpolation: d⋅b + (1-d)a
            % Bound d to [0,1]
            a = T.States_OP(indOP1,1)';                 % <= only dep. on X
            b = T.States_OP(indOP2,1)';                 % <= only dep. on X
            c = T.posBase(iT,1)';                       % <= only dep. on X
            
            d = ((b-a)'*(c-a))/((b-a)'*(b-a));
            d = min(max(d,0),1);
            
            r1 = 1-d;
            r2 = d;
            
            intOPs{iT}(iiT,:) = [indOP1,r1,indOP2,r2];
        end
    end % every influencing turbine
end % every turbine

end

