function paramFLORIS = parameterFLORIS()
%PARAMETERFLORIS returns all parameters needed to evaluate the FLORIS model
Bart            = false;
param_clc_2024  = false;
if Bart
    paramFLORIS.alpha   = 1.088; %#ok<UNRCH>
    paramFLORIS.beta    = 0.222;
    paramFLORIS.k_a     = 0.537;
    paramFLORIS.k_b     = -0.000848;
    paramFLORIS.k_fa    = 7.84;
    paramFLORIS.k_fb    = 4.57;
    paramFLORIS.k_fc    = 0.43;
    paramFLORIS.k_fd    = -0.246;
    paramFLORIS.eta     = 1.0
    paramFLORIS.p_p     = 2.2;
    paramFLORIS.airDen  = 1.225;    % Air density kg/m^3 (SOWFA)
    paramFLORIS.TIexp   = 3;        % TI crosswind distribution is modelled as 
                                    %   a factor of sigma_y and sigma_z. This
                                    %   is in contrast to previous literature
                                    %   which assumed a fixed box area and the
                                    %   turbines therein as source of higher TI
elseif param_clc_2024
    paramFLORIS.alpha   = 0.7515;
    paramFLORIS.beta    = 0.3483;
    paramFLORIS.k_a     = 0.7701;
    paramFLORIS.k_b     = 0.003678;
    paramFLORIS.k_fa    = 0.73;
    paramFLORIS.k_fb    = 0.0467;
    paramFLORIS.k_fc    = 0.2055;
    paramFLORIS.k_fd    = -0.32;
    % Note: SOWFA ADM power is too high, so eta > 1 is justified, otherwise the wind speed is too high. 
    paramFLORIS.eta     = 1.39;
    paramFLORIS.p_p     = 2.5433;
    paramFLORIS.airDen  = 1.225;    % Air density kg/m^3 (SOWFA)
    paramFLORIS.TIexp   = 1.5374;        % TI crosswind distribution is modelled as 
                                    %   a factor of sigma_y and sigma_z. This
                                    %   is in contrast to previous literature
                                    %   which assumed a fixed box area and the
                                    %   turbines therein as source of higher TI
else
    paramFLORIS.alpha   = 2.32;
    paramFLORIS.beta    = 0.154;
    paramFLORIS.k_a     = 0.38371;
    paramFLORIS.k_b     = 0.003678;
    paramFLORIS.k_fa    = 0.73;
    paramFLORIS.k_fb    = 0.8325;
    paramFLORIS.k_fc    = 0.0325;
    paramFLORIS.k_fd    = -0.32;
    paramFLORIS.eta     = 1.0;
    paramFLORIS.p_p     = 2.2;
    paramFLORIS.airDen  = 1.225;    % Air density kg/m^3 (SOWFA)
    paramFLORIS.TIexp   = 3;        % TI crosswind distribution is modelled as 
                                    %   a factor of sigma_y and sigma_z. This
                                    %   is in contrast to previous literature
                                    %   which assumed a fixed box area and the
                                    %   turbines therein as source of higher TI
end

end

