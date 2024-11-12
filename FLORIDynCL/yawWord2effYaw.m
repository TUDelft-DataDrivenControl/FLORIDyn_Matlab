function effYaw = yawWord2effYaw(yawWorld,windDir)
%YAWWORD2EFFYAW returns the effective yaw angle in radians
effYaw = deg2rad(yawWorld - windDir);
end

