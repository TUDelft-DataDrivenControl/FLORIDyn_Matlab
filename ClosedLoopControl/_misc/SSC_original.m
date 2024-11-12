% clear all; clc;
addpath(genpath('FLORISSE_M'))

% Setup zeroMQ server
zmqServer = zeromqObj('/home/bmdoekemeijer/OpenFOAM/bmdoekemeijer-2.4.0/jeromq/jeromq-0.4.4-SNAPSHOT.jar',3321,3600,true);

% Load the yaw setpoint LUT and set-up a simple function
nTurbs = 9;

% Initial control settings
yawAngleArrayOut  	= 255.0*ones(1,nTurbs); % Initial settings
torqueArrayOut 		= 0.0 * ones(1,nTurbs); % Not used -- placeholder
pitchAngleArrayOut 	= 0.0 * ones(1,nTurbs);

dataSend = setupZmqSignal(torqueArrayOut,yawAngleArrayOut,pitchAngleArrayOut);
%dataSend = setupZmqSignal(yawAngleArrayOut,torqueArrayOut,pitchAngleArrayOut);

firstRun = true;
            
disp(['Entering wind farm controller loop...']);
while 1
    % Receive information from SOWFA
    dataReceived = zmqServer.receive();
    currentTime  = dataReceived(1,1);
    measurementData = dataReceived(1,2:end);

    % Measurements: [genPower,rotSpeedF,azimuth,rotThrust,rotTorque,genTorque,nacYaw,bladePitch]
    %generatorPowerArray = measurementVector(1:8:end);
    generatorPowerArray = measurementData(1:8:end);
    rotorSpeedArray     = measurementData(2:8:end);
    azimuthAngleArray   = measurementData(3:8:end);
    rotorThrustArray    = measurementData(4:8:end);
    rotorTorqueArray    = measurementData(5:8:end);
    genTorqueArray      = measurementData(6:8:end);
    nacelleYawArray     = measurementData(7:8:end);
    bladePitchArray     = measurementData(8:8:end);

    if firstRun
        dt = rem(currentTime,10)
        firstRun = false;
    end
    
    % Update message string
    t = currentTime;
    

    if and(t >= 20600, t<=20900)  
		% Turbines stay orthogonal to uniform 60 deg wind dir change) 
		yawAngleArrayOut = ones(1,nTurbs)*interp1([20600,20900],[255,195],t);
    end
	    
    disp([datestr(rem(now,1)) '__    Synthesizing message string.']);
    dataSend = setupZmqSignal(torqueArrayOut,yawAngleArrayOut,pitchAngleArrayOut);
    % dataSend = setupZmqSignal(yawAngleArrayOut,torqueArrayOut,pitchAngleArrayOut); %torque went to 0
       
    % Send a message (control action) back to SOWFA
    zmqServer.send(dataSend);
end

% Close connection
zmqServer.disconnect()

function [dataOut] = setupZmqSignal(torqueSignals,yawAngles,pitchAngles)
	dataOut = [];
    for i = 1:length(yawAngles)
        dataOut = [dataOut torqueSignals(i) yawAngles(i) pitchAngles(i)];
    end
end