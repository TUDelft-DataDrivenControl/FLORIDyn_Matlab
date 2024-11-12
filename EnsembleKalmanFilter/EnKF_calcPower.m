function [P,EnKF] = EnKF_calcPower(C_u, EnKF, TD, paramFLORIS, TStartI)
%ENKF_CALCPOWER Calculates the power of the turbines
%   This is done as the wind field states have been projected on a common
%   OP location, at which we correct the states.

% Free wind speed from the ensembles
Vel = EnKF.States.Vel;
nT  = length(TD);
red = zeros(nT, EnKF.nE);
yaw = deg2rad(EnKF.States_T(TStartI,2:EnKF.nStatesT:end));

for iE = 1:EnKF.nE
    red(:,iE) = EnKF.M{iE}.("Foreign Reduction [%]")(end-nT+1:end)...
        * 10^(-2);
end

ueff = red .* (C_u * Vel);
a    = 1/3;
Cp   = 4*a.*(1-a).^2;
P    = 0.5* paramFLORIS.airDen*(TD/2).^2....
        .* pi .* Cp .* ueff.^3 .* paramFLORIS.eta * 10^(-6) .* ...
        cos(yaw).^paramFLORIS.p_p ;

% Relace the measurements for the new ones
EnKF = EnKF_overwritePowerM(EnKF,nT,P);

end

