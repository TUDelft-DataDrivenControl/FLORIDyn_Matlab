function [] = animateTWF()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Parameters and constants
nOP = 30;
dt  = 10;
u   = 8;
D   = 175; 
T1  = [0,0]';
T2  = [5*D,0]';
OPs_0 = repmat(T1',nOP,1) + (1:30)'; %[x,y]
OPs_1 = [(1:30)',(1:30)']; %[x,y]
i_Change1 = 50;
i_Change2 = 100;
%% Function handles
% Surrogate deflection function (made up for demonstration)
delta =@(x1) tanh(0.01*x1)+sqrt(0.02*x1);
% Arbitrary wind direction change
dir =@(t) 10 + 0.5*(tanh(0.3*(t-i_Change1))+1)*30 ...
    - 0.5*(tanh(0.3*(t-i_Change2))+1)*30;

% Eq.(2) Rotational Matrix 
R_01 =@(phi) [cos(phi),-sin(phi);sin(phi),cos(phi)];

f   = figure(1000);
%% Simulation
for i = 1:200
    % Get new wind direction
    phi = dir(i);
    
    % OP shift
    OPs_0(2:end,:)  = OPs_0(1:end-1,:);
    OPs_1(2:end,:)  = OPs_1(1:end-1,:);
    OPs_0(1,:)      = T1';
    OPs_1(1,:)      = 0;
    % Store old K_1 coordinates
    OPs_1_k = OPs_1;
    
    % Eq.(3)
    OPs_1(:,1) = OPs_1(:,1) + dt*u;
    % Eq.(4)
    OPs_1(:,2) = delta(OPs_1(:,1));
    % Eq.(1)
    OPs_0 = (OPs_0' + R_01(deg2rad(phi))*(OPs_1'-OPs_1_k'))';
    
    % Set Location, here turbine 2
    L_0 = T2;
    
    % Eq.(5) and (6)
    [OP_star_0,OP_star_1] = findOPstar(OPs_0,OPs_1,L_0);
    
    % Eq.(7)
    L_T1_1 = OP_star_1 + R_01(deg2rad(phi))'*(L_0 - OP_star_0);
    
    % Eq.(8) in such a way that it fits the wind dir of K_0
    t_T1_2 = L_0 - R_01(deg2rad(phi))*L_T1_1;
    
    %% Plot
    plot(OPs_0(:,1),OPs_0(:,2),'--o')
    hold on
    scatter(OP_star_0(1),OP_star_0(2),'filled')
    scatter(T2(1),T2(2),'filled')
    scatter(t_T1_2(1),t_T1_2(2),'filled')
    scatter(T1(1),T1(2),'filled')
    quiver(800,-100,100*cos(deg2rad(phi)),100*sin(deg2rad(phi)),0,...
        'k','LineWidth',2,'AutoScale','off',...
            'MaxHeadSize',100)
    plot([T2(1),OP_star_0(1),t_T1_2(1),T2(1)],...
        [T2(2),OP_star_0(2),t_T1_2(2),T2(2)],':')
    hold off
    axis equal
    legend('OPs','OP^*','T1','T1_2','T2','Wind dir','Location','northwest')
    grid on
    title(['t=' num2str(i*dt) ' s, wind dir. change around '...
        num2str(i_Change1*dt) ' s and ' num2str(i_Change2*dt) ' s.'])
    xlim([-100,1000])
    ylim([-300,800])
    drawnow
    f.Position(3) = 496;
    f.Position(4) = 496;
    saveas(f,['TWFAnimation/' num2str(i) '.png'])
end


end

function [OP_star_0,OP_star_1] = findOPstar(OPs_0,OPs_1,L_0)
nOP = size(OPs_0,1);
% Get distance
dist = sqrt(sum((OPs_0 - L_0').^2,2));

% Find index of min
[~,I] = min(dist);

if and(I<nOP,I>1)
    % there is a previous and following OP
    if dist(I+1)>dist(I-1)
        % previous is closer
        x_OP2_0 = OPs_0(I-1,:)';
        x_OP2_1 = OPs_1(I-1,:)';
    else
        % following is closer
        x_OP2_0 = OPs_0(I+1,:)';
        x_OP2_1 = OPs_1(I+1,:)';
    end
elseif I==nOP
    % Last OP
    x_OP2_0 = OPs_0(I-1,:)';
    x_OP2_1 = OPs_1(I-1,:)';
else
    % First OP
    x_OP2_0 = OPs_0(I+1,:)';
    x_OP2_1 = OPs_1(I+1,:)';
end
x_OP1_0 = OPs_0(I,:)';
x_OP1_1 = OPs_1(I,:)';

% Eq. (5) 
w = ((x_OP2_0 - x_OP1_0)'*(L_0-x_OP1_0))/...
    ((x_OP2_0 - x_OP1_0)'*(x_OP2_0 - x_OP1_0));

% Limit for first and last OP
w = min(1,max(0,w));

% Eq.(6) 
OP_star_0 = (1-w)*x_OP1_0 + w*x_OP2_0;
OP_star_1 = (1-w)*x_OP1_1 + w*x_OP2_1;
end

