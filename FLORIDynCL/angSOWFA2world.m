function rad_World = angSOWFA2world(deg_SOWFA)
%ANGSOWFA2WORLD 
% Angle conversion SOWFA to world coordinates
% deg_F = -deg_S + 270 deg
% yaw angle defined clockwise, but for calculations counterclockwise
deg_World = 270 - deg_SOWFA;
rad_World = deg2rad(deg_World);
end

