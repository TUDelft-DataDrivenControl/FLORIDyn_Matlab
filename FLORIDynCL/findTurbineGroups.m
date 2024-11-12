function Tdep = findTurbineGroups(T,paramFLORIDyn)
%FINDTURBINEGROUPS determines which turbines are influencing the i-th
%turbine. All turbines +10D upstream and +-3D crossstream are considered. 
% =============== INPUTS ================================================ %
%   T       Turbine OP and wind field data
%   paramFLORIDyn
%       .deltaDW  Down Wind Diameters distance from closest OP to be influenced
%       .deltaCW  Cross Wind Diameters distance from closest OP t.b.infl.
%       .deltaUW  Up Wind Diameters distance from closest OP t.b.infl.
% =============== OUTPUTS =============================================== %
%   Tdep    Dependencies
% ======================================================================= %
dw = paramFLORIDyn.deltaDW;
cw = paramFLORIDyn.deltaCW;
uw = paramFLORIDyn.deltaUW;


Tdep = cell(T.nT,1);
dep = false(T.nT);
R01 =@(phi) [cos(phi) sin(phi); -sin(phi) cos(phi)];
for iT = 1:T.nT
    for iiT = 1:T.nT
        % for every other turbine find out if they less than the threshold
        % away from the centerline
        if iiT == iT; continue; end
        
        % find closest OP from the wake of turbine iT to turbine iiT and 
        %   extract index
        distOP_iiT = sum((...
            T.posBase(iiT,1:2) - T.States_OP((0:T.nOP-1)+T.StartI(iT),1:2)...
            ).^2,2);
        [~,I_op] = min(distOP_iiT);
        
        % Calculate the vector from the Turbine to the OP
        I_op = I_op + T.StartI(iT) - 1;
        phi  = angSOWFA2world(T.States_WF(I_op,2));
        r0   = T.States_OP(I_op,1:2) - T.posBase(iiT,1:2);

        % Roate the vector in down-/crossstream coordinates
        r1   = R01(phi)*r0';
        
        % r1(1) downstream distance: pos -> OP is downstream
        %                            neg -> OP is upstream
        % r1(2) crosstream distance
        if and(and( ...                     % Distance checks
                -r1(1) <= uw*T.D(iT), ...   % -> Upstream distance 
                 r1(1) <= dw*T.D(iT)), ...  % -> Downstream distance
                abs(r1(2)) <= cw*T.D(iT)... % -> Crosswind distance
                )
            dep(iiT,iT) = true;
        end

        %% Old caluclation
        % r1 = abs(R01(phi)*r0');
        
        % Note for future
        %   In this code the up-/down-stream distance and cross-stream
        %   distance are checkt to determine if there might be an impact of
        %   the OP of turbine iT onto turbine iiT. r(1) is the up-/down-
        %   stream distance, r(2) cross-stream. Since the absolute is
        %   taken, and the distance is chosen generously, dependencies can
        %   be suggested where a downstream wake would impact an upstream
        %   turbine. This is prevented at a later stage, but the fact that
        %   it is recorded can cause issues if the 'influence' correction
        %   is used for wind field variables.

        % if and(r1(1)<=dw*T.D(iT),r1(2)<=cw*T.D(iT)) %and(and(abs(r1(1))<=dw*T.D(iT), r1(1)<0), abs(r1(2))<=cw*T.D(iT)) 
        %     dep(iiT,iT) = true;
        % end
    end
end

for iT = 1:T.nT
    Tdep{iT} = find(dep(iT,:));
end
end

