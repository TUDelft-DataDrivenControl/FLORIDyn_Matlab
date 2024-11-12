function states = States()
%STATES States of the FLORIS model
%   T  -> turbine states needed to calculate the wake
%   OP -> Observation point states (world and wake coordinates)
%   WF -> Wind field states
% ======================================================================= %
states.T_names  = {'a','yaw','TI'};
states.Turbine  = length(states.T_names);

states.OP_names = {'x0','y0','z0','x1','y1','z1'};
states.OP       = length(states.OP_names);

states.WF_names = {'wind_vel','wind_dir','TI0'};
states.WF       = length(states.WF_names);
end

