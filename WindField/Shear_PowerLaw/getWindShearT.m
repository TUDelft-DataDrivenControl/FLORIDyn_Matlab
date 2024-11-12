function shear = getWindShearT(WindShear,z_norm)
%GETWINDSHEART Return the shear factor u_eff = shear * u_referenceHight
% POWER LAW implementation
%   expects a WindShearPowerLaw.csv with the shear coefficient
% 
% ======================================================================= %
% WindShear = Holds shear coefficient and reference height
%          .z0      = reference height
%          .alpha   = shear coefficient
% z         = height(s)
% ======================================================================= %
shear = (z_norm).^WindShear.alpha;
end

