function Con = controller(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS,CLC)
%CONTROLLER determines the control inputs for the near future and beyond.
%   It utilizes a given cost-function to evaluate the 
%% Brute force Bezier_yaw
% This controller is based on Bezier-Curves, which allow a smooth
% transition between distinct setpoints. These setpoints are subject to the
% optimization.
% Bezier-curves further allow the user to calculate the derivative over
% time, which allows us to enforce a yaw-rate limit (or axial induction)

%% The brute force part
%   This controller essentially performs a grid-search for the best option
%   and maps out the solution space. It is meant to learn about the cost
%   function and the nature of the optimization problem.


end

