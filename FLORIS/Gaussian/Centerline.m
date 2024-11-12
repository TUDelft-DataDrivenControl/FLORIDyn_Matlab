function delta = Centerline(States_OP,States_T,States_WF,paramFLORIS,D)
%CENTERLINE model

%% Parameters
k_a     = paramFLORIS.k_a;
k_b     = paramFLORIS.k_b;
alpha   = paramFLORIS.alpha;
beta    = paramFLORIS.beta;

%% States
C_T     = CalcCt(States_T(:,1),States_T(:,2));
yaw     = -deg2rad(States_T(:,2));
I       = sqrt(States_T(:,3).^2+States_WF(:,3).^2); % I_f & I_0
OPdw    = States_OP(:,4);

%% Calc x_0 (Core length)
% [1] Eq. 7.3
x_0 = (cos(yaw).*(1+sqrt(1-C_T))./...
    (sqrt(2)*(alpha*I+beta*(1-sqrt(1-C_T))))).*D;

%% Calc k_z and k_y based on I
%[2] Eq.8
k_y = k_a*I + k_b;
k_z = k_y;

%% Get field width y
% [1] Eq. 7.2
% To fit the field width, the value linearly increases from 0 to max for dw
% positions before x_0
zs = zeros(size(OPdw));
sig_y = ...
    max([OPdw-x_0,zs],[],2)   .* k_y +...
    min([OPdw./x_0,zs+1],[],2).* cos(yaw) .* D/sqrt(8);

%% Get field width z
% [1] Eq. 7.2
sig_z = ...
    max([OPdw-x_0,zs],[],2)    .* k_z +...
    min([OPdw./x_0,zs+1],[],2) .* D/sqrt(8);

%% Calc Theta
%[1] Eq. 6.12
Theta = 0.3*yaw./cos(yaw).*(1-sqrt(1-C_T.*cos(yaw)));

%% Calc Delta / Deflection
%[1] Eq. 7.4 (multiplied with D and the second part disabled for dw<x_0)
%   Part 1 covering the linear near field and constant for the far field
delta_nfw  = Theta.*min([OPdw,x_0],[],2);

%   Part 2, smooth angle for the far field, disabled for the near field.
delta_fw_1 = Theta/14.7.*sqrt(cos(yaw)./(k_y.*k_z.*C_T)).*(2.9+1.3*sqrt(1-C_T)-C_T);
delta_fw_2 = log(...
    (1.6+sqrt(C_T)).*...
    (1.6.*sqrt((8*sig_y.*sig_z)./(D.^2.*cos(yaw)))-sqrt(C_T))./(...
    (1.6-sqrt(C_T)).*...
    (1.6.*sqrt((8*sig_y.*sig_z)./(D.^2.*cos(yaw)))+sqrt(C_T))...
    ));

% Combine deflection parts
deltaY = delta_nfw + ...
    (sign(OPdw-x_0)/2+0.5).*...
    delta_fw_1.*delta_fw_2.*D;


% Deflection in y and z direction
delta = [deltaY, zeros(size(deltaY))];
end

