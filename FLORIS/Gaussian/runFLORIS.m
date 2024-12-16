function [T_red_arr,T_aTI_arr, T_Ueff, T_weight] = runFLORIS(LocationT, States_WF, States_T, D, paramFLORIS,WindShear)
%RUNFLORIS Receives the interpolated info from FLORIDyn and returns the
% measurements for the last turbine
% OUTPUT
%   T_red_arr   := Reduction by every turbine and environment (last entry)
%   T_aTI_arr   := Added turbulence intensity by the turbines
%   T_Ueff      := Effective wind speeds

% Get rotor points of interest
if D(end)>0
    [RPl,RPw] = discretizeRotor(paramFLORIS.RotorPoints);
else
    % Plotting grid point
    RPl = [0,0,0];
    RPw = 1;
end

% Move them according to yaw and nacelle location of the last turbine
tmp_yaw = deg2rad(States_T(end,2));
% Rotation matrix based on yaw angle
R       = [cos(tmp_yaw),sin(tmp_yaw),0; -sin(tmp_yaw),cos(tmp_yaw),0;0,0,1];
% Scaled rotor point locations
RPl     = (R*(RPl.*D(end))')'+ LocationT(end,:);

if length(D) == 1
    % Single turbine wind farm -> only wind shear influence
    redShear    = getWindShearT(WindShear,RPl(:,3)/LocationT(3));
    T_red_arr   = RPw'*redShear;
    T_aTI_arr   = [];
    T_Ueff      = [];
    return
end
    
% Preallocate outputs
T_red_arr = ones(length(D),1); % nT-1 turbines and one environment entry
T_aTI_arr = zeros((length(D)-1),1);
T_weight  = zeros((length(D)-1),1);

% Calculate impact by each turbine
for iT = 1:(length(D)-1)
    % Retrieve the orientation
    if size(States_WF,2) == 4
        % Original orientation
        tmp_phi = angSOWFA2world(States_WF(iT,4));
    else
        % Wind dir
        tmp_phi = angSOWFA2world(States_WF(iT,2));
    end
    
    % Transform rotor points to the wake coordinate system of iT
    tmp_RPs = RPl - LocationT(iT,:);
    R = [cos(tmp_phi),sin(tmp_phi),0;-sin(tmp_phi),cos(tmp_phi),0;0,0,1];
    tmp_RPs = (R*tmp_RPs')';
    
    if tmp_RPs(1)<=10
        % influenced Turbine is upstream -> no influence
        continue;
    end
    % =============== DEBUG PLOT ======================
%     figure
%     subplot(2,1,1)
%     scatter3(LocationT(end,1),LocationT(end,2),LocationT(end,3),'filled')
%     hold on
%     scatter3(RPl(:,1),RPl(:,2),RPl(:,3),'filled')
%     scatter3(LocationT(iT,1),LocationT(iT,2),LocationT(iT,3))
%     hold off
%     axis equal
%     grid on
%     title('World coordinate relations')
%     subplot(2,1,2)
%     scatter3(tmp_RPs(:,1),tmp_RPs(:,2),tmp_RPs(:,3),'filled')
%     hold on
%     scatter3(0,0,0)
%     axis equal
%     hold off
%     title('Wake coordinate relations')
    
    % Turbine states (iT)
    a   = States_T(iT,1);
    yaw = -deg2rad(States_T(iT,2));
    TI  = States_T(iT,3);
    Ct  = CalcCt(a,States_T(iT,2));
    TI0 = States_WF(iT,3);
    
    % Calculate variables for Gaussian model
    [sig_y, sig_z, C_T, x_0, delta, pc_y, pc_z] = ...
        getVars(tmp_RPs,a,Ct,yaw,TI,TI0,paramFLORIS,D(iT));
    % create an radius value of the core and cw values and figure out if the
    % OPs are in the core or not
    cw_y    = tmp_RPs(:,2) - delta(:,1);
    cw_z    = tmp_RPs(:,3) - delta(:,2);
    
    phi_cw  = atan2(cw_z,cw_y);
    r_cw    = sqrt(cw_y.^2 + cw_z.^2);
    core    = or(...
        r_cw < sqrt((0.5*pc_y.*cos(phi_cw)).^2 + ...
        (0.5*pc_z.*sin(phi_cw)).^2), ...
        tmp_RPs(:,1)==0);
    nw = tmp_RPs(:,1)<x_0;
    
    % Calculate reduction r
    tmp_RPs_r = zeros(size(RPw));
    tmp_RPs_r(core) = 1-sqrt(1-C_T);
    
    % Remove core from crosswind pos and calculate speed reduction
    fw = ~nw;
    gaussAbs = zeros(size(core));
    
    gaussAbs(nw) = 1-sqrt(1-C_T);
    gaussAbs(fw) = 1-sqrt(1-C_T...
        .*cos(yaw)./(8*(sig_y(fw).*sig_z(fw)./D(iT).^2)));
    
    gaussWght = ones(size(core));
    gaussWght(~core) = ...
        exp(-0.5.*((cw_y(~core)-cos(phi_cw(~core)).*pc_y(~core)*0.5)./sig_y(~core)).^2).*...
        exp(-0.5.*((cw_z(~core)-sin(phi_cw(~core)).*pc_z(~core)*0.5)./sig_z(~core)).^2);
    tmp_RPs_r(~core) = gaussAbs(~core).*gaussWght(~core);

    T_weight(iT) = sum(gaussWght);

    % Combine reductions
    T_red_arr(iT) = (1-RPw'*(tmp_RPs_r));
    
    % Added turbulence levels (squares summed)
    %   Calculate influence
    T_addedTI_tmp = (...
        paramFLORIS.k_fa*(...
        a.^paramFLORIS.k_fb .* ...
        TI0.^paramFLORIS.k_fc .* ...
        (mean(tmp_RPs(:,1))./D(iT)).^paramFLORIS.k_fd)...
        );
    %   Weight influence
    %   cw_y: y crosswind distance of the RPs to the deflected centreline
    %   pc_y: y crosswind width of the potential core (Only relevant in
    %           near field)
    T_aTI_arr(iT) = T_addedTI_tmp * RPw'* ...
        (exp(-0.5.*((cw_y-cos(phi_cw).*pc_y*0.5)./...
        (paramFLORIS.TIexp*sig_y)).^2).*...
        exp(-0.5.*((cw_z-sin(phi_cw).*pc_z*0.5)./...
        (paramFLORIS.TIexp*sig_z)).^2));
    
end
redShear = getWindShearT(WindShear,RPl(:,3)/LocationT(iT,3));
T_red_arr(end) = RPw'*redShear;

T_red       = prod(T_red_arr);
T_Ueff      = States_WF(end,1)*T_red;
end

