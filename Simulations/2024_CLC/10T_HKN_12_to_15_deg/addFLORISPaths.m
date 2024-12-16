%% Links to the specific version of FLORIS which is used
addpath(['.' filesep 'FLORIS' filesep Sim.FLORIS]);
addpath(['.' filesep 'WindField' filesep 'Velocity_' Wind.Input_Vel]);
addpath(['.' filesep 'WindField' filesep 'Direction_' Wind.Input_Dir]);
addpath(['.' filesep 'WindField' filesep 'TI_' Wind.Input_TI]);
addpath(['.' filesep 'WindField' filesep 'Shear_' Wind.Input_Shr]);
addpath(['.' filesep 'FLORIS' filesep 'Discretization' filesep Sim.RotorDiscret]);
addpath(['.' filesep 'Controller' filesep 'Yaw_' Con.Yaw]);
addpath(['.' filesep 'Correction' filesep 'Direction_' Wind.Correction.Dir]);
addpath(['.' filesep 'Correction' filesep 'Velocity_' Wind.Correction.Vel]);
addpath(['.' filesep 'Correction' filesep 'TI_' Wind.Correction.TI]);
p = mfilename('fullpath');
Sim.PathToSim = p(1:end-14);
clear p