function P = getPower(T,M,paramFLORIS,Con)
a   = T.States_T(T.StartI,1);
%yaw = yawWord2effYaw(T.States_T(T.StartI,2),T.States_WF(T.StartI,2));
yaw = deg2rad(T.States_T(T.StartI,2));

Cp      = 4*a.*(1-a).^2;
ueff    = M(:,3);

if Con.tanhYaw
    P = 0.5*paramFLORIS.airDen*(T.D/2).^2....
        *pi.*Cp.*ueff.^3.* paramFLORIS.eta.*...
        cos(yaw).^paramFLORIS.p_p.*...
        (0.5*tanh((-yaw+deg2rad(Con.yawRangeMax))*50)+.5) .* ...
        (-0.5*tanh((-yaw+deg2rad(Con.yawRangeMin))*50)+.5);
else
    P = 0.5*paramFLORIS.airDen*(T.D/2).^2....
        *pi.*Cp.*ueff.^3.* paramFLORIS.eta.*...
        cos(yaw).^paramFLORIS.p_p;
end
end