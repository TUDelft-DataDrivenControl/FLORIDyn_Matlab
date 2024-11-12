function Ct = CalcCt(a,~)
%CALCCT Calculates Ct based on the axial induction factor and yaw angle
% INPUT
%   a   : Axial indcution factor
%   yaw : Yaw angle (deg, counter clockwise)
Ct = 4 * a .* (1 - a);
end

