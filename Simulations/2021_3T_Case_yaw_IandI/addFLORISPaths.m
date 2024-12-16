%% Links to the specific version of FLORIS which is used
% NOTE: There are no values to be modified in this file, it has to be in
% the simulation folder to automatically retreive the path to the
% simulation folder and store it in "Sim.PathToSim".
%
% ======= DO NOT MODIFY ================================================= %
if ispc
    % Microsoft pathing
    addpath(['.\FLORIS\' Sim.FLORIS]);
    addpath(['.\WindField\Velocity_' Wind.Input_Vel]);
    addpath(['.\WindField\Direction_' Wind.Input_Dir]);
    addpath(['.\WindField\TI_' Wind.Input_TI]);
    addpath(['.\WindField\Shear_' Wind.Input_Shr]);
    addpath(['.\FLORIS\Discretization\' Sim.RotorDiscret]);
    addpath(['.\Controller\Yaw_' Con.Yaw]);
    addpath(['.\Correction\Direction_' Wind.Correction.Dir]);
    addpath(['.\Correction\Velocity_' Wind.Correction.Vel]);
    addpath(['.\Correction\TI_' Wind.Correction.TI]);
else
    % Unix pathing
    addpath(['FLORIS/' Sim.FLORIS]);
    addpath(['WindField/Velocity_' Wind.Input_Vel]);
    addpath(['WindField/Direction_' Wind.Input_Dir]);
    addpath(['WindField/TI_' Wind.Input_TI]);
    addpath(['WindField/Shear_' Wind.Input_Shr]);
    addpath(['FLORIS/Discretization/' Sim.RotorDiscret]);
    addpath(['Controller/Yaw_' Con.Yaw]);
    addpath(['Correction/Direction_' Wind.Correction.Dir]);
    addpath(['Correction/Velocity_' Wind.Correction.Vel]);
    addpath(['Correction/TI_' Wind.Correction.TI]);
end
p = mfilename('fullpath');
Sim.PathToSim = p(1:end-14);
clear p
% ======= DO NOT MODIFY ================================================= %