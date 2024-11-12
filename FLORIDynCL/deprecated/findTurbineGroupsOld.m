function Tdep = findTurbineGroupsOld(Tpos, D, nT, WF, uw, cw)
%FINDTURBINEGROUPS determines which turbines are influencing the i-th
%turbine. All turbines +10D upstream and +-3D crossstream are considered. 
% =============== INPUTS ================================================ %
%   Tpos    Turbine world coordinates only xy used
%   D       Turbine Diameters
%   nT      Number of turbines
%   WF      Wind field states at the turbine locations
%   uw      Up wind Diameters distance
%   cw      Crowss wind Diameters distance
% =============== OUTPUTS =============================================== %
%   Tdep    Dependencies
% ======================================================================= %
%           T.dep = findTurbineGroups(T.posBase, T.D, T.nT, ...
%                 T.States_WF(T.StartI,:), paramFLORIDyn.deltaDW,...
%                 paramFLORIDyn.deltaCW);
% ======================================================================= %
Tdep = cell(nT,1);
for iT = 1:nT
    % Get wind direction at turbine location
    phi = angSOWFA2world(WF(iT,2));
    
    % Rotate wind farm accordingly and remove offset from current turbine
    xT = cos(phi).*Tpos(:,1) + sin(phi).*Tpos(:,2);
    yT = -sin(phi).*Tpos(:,1) + cos(phi).*Tpos(:,2);
    xT = xT - xT(iT);
    yT = abs(yT - yT(iT));
    
    % Determine indeces of turbines 
    % ///////// TMP Disable -> causes sharp cuts, has to be rethought
    % Better approach than turbine searches influeces is that influences
    % inform turbine (Distance from closest OP smaller than ..)
%     dep = and(...
%             and(xT<0,xT>-uw*D(iT)),...
%             yT<=cw*D(iT));
    dep = true(size(xT));
    dep(iT) = false;
    Tdep{iT} = find(dep);
end
end

