function paramFLORIS = parameterFLORIS()
%PARAMETERFLORIS returns all parameters needed to evaluate the FLORIS model
Bart = false;
if Bart
    paramFLORIS.alpha   = 1.088; %#ok<UNRCH>
    paramFLORIS.beta    = 0.222;
    paramFLORIS.k_a     = 0.537;
    paramFLORIS.k_b     = -0.000848;
    paramFLORIS.k_fa    = 7.84;
    paramFLORIS.k_fb    = 4.57;
    paramFLORIS.k_fc    = 0.43;
    paramFLORIS.k_fd    = -0.246;
else
    paramFLORIS.alpha   = 2.32;
    paramFLORIS.beta    = 0.154;
    paramFLORIS.k_a     = 0.38371;
    paramFLORIS.k_b     = 0.003678;
    paramFLORIS.k_fa    = 0.73;
    paramFLORIS.k_fb    = 0.8325;
    paramFLORIS.k_fc    = 0.0325;
    paramFLORIS.k_fd    = -0.32;
end

paramFLORIS.eta     = 1;            %0.7918; %0.8572;
paramFLORIS.p_p     = 2.2;
paramFLORIS.airDen  = 1.225; % Air density kg/m^3 (SOWFA)
end

