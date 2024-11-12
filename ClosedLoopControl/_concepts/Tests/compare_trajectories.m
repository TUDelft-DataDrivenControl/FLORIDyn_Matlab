nHorizonSteps = 50;
nControlSteps = 50;
nUpdateSteps = 20;
nIterations = 10;
nBranches = 1;
deltaT = 4;
rl = .3;

rn = rl*nControlSteps*deltaT;
tn = linspace(0,1,nControlSteps+1);
%% Part 1 - bezier approach
%  Cubic bezier curves
%   continuous transition

% Starting conditions
gamma_d = 0;
gamma_0 = 0;
t_0 = 0;

% Bezier trajectory for normalised time
bc = @(t_n,g0,g1,g2,g3) ...
    (1-t_n).^3*g0 +...
    3*(1-t_n).^2 .* t_n * g1 +...
    3*(1-t_n).* t_n.^2 * g2 +...
    t_n.^3 * g3;

% Derivative of the Bezier trajectory for normalised time
bc_d = @(t_n,g0,g1,g2,g3) ...
    3*(1-t_n).^2  * (g1 - g0)  +...
    6*(1-t_n).* t_n * (g2 - g1) +...
    3*t_n.^2 * (g3 - g2);

% Translation of a opimization value [0,1] into g2
o1_to_g2 = @(o1, g1, rn) (1-o1)*(...
    1/3*( -sqrt( rn^2 + ...
    3*rn * g1)-...
    rn + 3 * g1)) + ...
    o1 * (...
    -1/3*( -sqrt( rn^2 - ...
    3*rn * g1)-...
    rn - 3 * g1));

%% Part 1.1 Bezier curve approach with 1 dof
figure
subplot(2,2,1)
hold on



for nI = 1:nIterations
    % generate n seeds
    o1 = rand(nBranches,1);

    for iO = 1:nBranches
        % calculate trajectory
        g2 = o1_to_g2(o1(iO), gamma_d/3, rn);
        tr = bc(tn,0,gamma_d/3,g2,g2) + gamma_0;

        if iO == nBranches
            % Plot chosen one
            plot((0:nControlSteps)*deltaT+t_0, tr,'--',"Color",...
                [0 184 200]./255, "LineWidth",1.5)
            plot((nControlSteps:nHorizonSteps)*deltaT+t_0,...
                tr(end)*ones(size(nControlSteps:nHorizonSteps)),'--',"Color",...
                [0 184 200]./255, "LineWidth",1.5)
            plot((0:nUpdateSteps-1)*deltaT+t_0, tr(1:nUpdateSteps),"Color",...
                [12 35 64]./255, "LineWidth",2.5)

            % Add line to mark end of the control horizon.
            plot([nControlSteps, nControlSteps]*deltaT + t_0,...
                [-2.9,2.9]+tr(nControlSteps),...
                "LineWidth",2,"Color",[255,184,28]./255)

            % Assign new starting conditions
            gamma_0 = bc(tn(nUpdateSteps),0,gamma_d/3,g2,g2) + gamma_0;
            gamma_d = bc_d(tn(nUpdateSteps),0,gamma_d/3,g2,g2);
            t_0 = t_0 + (nUpdateSteps-1)*deltaT;
        else
            % Plot disregarded ones
            plot((0:nControlSteps)*deltaT+t_0,tr,"Color", ...
                [194 194 194]./255, "LineWidth",1)
            plot((nControlSteps:nHorizonSteps)*deltaT+t_0,...
                tr(end)*ones(size(nControlSteps:nHorizonSteps)),"Color",...
                [194 194 194]./255, "LineWidth",1)
        end
    end
    
    % Add line to mark transition between two trajectories.
    plot([t_0, t_0], [-2.9,2.9]+gamma_0,...
        "LineWidth",2,"Color",[237,104,66]./255)
end
hold off
xlabel('time (s)')
ylabel('Orientation (deg)')
xlim([0,450])

clear bc bc_d g2 gamma_0 gamma_d iO nI o1_to_g2 

%% Bezier curve approach with 1 dof & time offset
% Starting conditions
gamma_d = 0;
gamma_0 = 0;
t_0 = 0;

% Bezier trajectory for normalised time
bc = @(t_n,g0,g1,g2,g3) ...
    (1-t_n).^3*g0 +...
    3*(1-t_n).^2 .* t_n * g1 +...
    3*(1-t_n).* t_n.^2 * g2 +...
    t_n.^3 * g3;

% Derivative of the Bezier trajectory for normalised time
bc_d = @(t_n,g0,g1,g2,g3) ...
    3*(1-t_n).^2  * (g1 - g0)  +...
    6*(1-t_n).* t_n * (g2 - g1) +...
    3*t_n.^2 * (g3 - g2);

% Translation of a opimization value [0,1] into g2
o1_to_g2 = @(o1, g1, rn) (1-o1)*(...
    1/3*( -sqrt( rn^2 + ...
    3*rn * g1)-...
    rn + 3 * g1)) + ...
    o1 * (...
    -1/3*( -sqrt( rn^2 - ...
    3*rn * g1)-...
    rn - 3 * g1));

% Calculation of the minimum time
%t_cn = @(g_d,r_n) min(abs(g_d)/r_n,1);

% Calculation of the start time
t_sn = @(g_d,r_n,o_t) o_t * (1 - min(abs(g_d)/r_n,1));

% new rate limit, based on time window
r_nstar = @(r_n,t_cn) r_n * t_cn;

% new angle based on time window
g_dstar = @(g_dn,t_cn) g_dn * t_cn; 


%% Plot BC w 1 dof & time offset
subplot(2,2,2)
hold on

for nI = 1:nIterations
    % generate n seeds
    o1 = rand(nBranches,1); % t_so
    o2 = rand(nBranches,1); % g_so

    for iO = 1:nBranches
        % Calculate start time & rescale rate limit and angle
        tsn = t_sn(gamma_d,rn,o1(iO));
        % Scale tsn such that it aligns with resolution of time steps
        tsn = tsn - mod(tsn,1/nControlSteps);
        r_ns = r_nstar(rn,1-tsn);
        g_ds = g_dstar(gamma_d,1-tsn);

        % calculate BC trajectory in control window
        g2 = o1_to_g2(o2(iO), g_ds/3, r_ns);
        tr = zeros(size(tn));
        tr(tn<tsn)  = tn(tn<tsn)*gamma_d + gamma_0;
        tr(tn>=tsn) = bc((tn(tn>=tsn)-tsn)/(1-tsn),...
            0,g_ds/3,g2,g2) + tsn*gamma_d + gamma_0;

        if iO == nBranches
            % Plot chosen one
            plot((0:nControlSteps)*deltaT+t_0, tr,'--',"Color",...
                [12 35 64]./255, "LineWidth",1.5)
            plot((nControlSteps:nHorizonSteps)*deltaT+t_0,...
                tr(end)*ones(size(nControlSteps:nHorizonSteps)),'--',"Color",...
                [12 35 64]./255, "LineWidth",1.5)
            plot((0:nUpdateSteps-1)*deltaT+t_0, tr(1:nUpdateSteps),"Color",...
                [0 184 200]./255, "LineWidth",2.5)

            % Add line to mark end of the control horizon.
            plot([nControlSteps, nControlSteps]*deltaT + t_0,...
                [-2.9,2.9]+tr(nControlSteps),...
                "LineWidth",2,"Color",[255,184,28]./255)

            % Assign new starting conditions
            gamma_0 = tr(nUpdateSteps);
            if nUpdateSteps/nControlSteps>tsn
                gamma_d = bc_d((nUpdateSteps/nControlSteps-tsn)/(1-tsn),...
                            0,g_ds/3,g2,g2)/(1-tsn);
            end
            %gamma_d = bc_d(tn(nUpdateSteps),0,gamma_d/3,g2,g2); % TODO fix
            t_0 = t_0 + (nUpdateSteps-1)*deltaT;
        else
            % Plot disregarded ones
            plot((0:nControlSteps)*deltaT+t_0,tr,"Color", ...
                [194 194 194]./255, "LineWidth",1)
            plot((nControlSteps:nHorizonSteps)*deltaT+t_0,...
                tr(end)*ones(size(nControlSteps:nHorizonSteps)),"Color",...
                [194 194 194]./255, "LineWidth",1)
        end
    end
    
    % Add line to mark transition between two trajectories.
    plot([t_0, t_0], [-2.9,2.9]+gamma_0,...
        "LineWidth",2,"Color",[237,104,66]./255)
end
hold off
xlabel('time (s)')
ylabel('Orientation (deg)')
xlim([0,450])

%clear bc bc_d g2 gamma_0 gamma_d iO nI o1_to_g2 

%% Part 2 Bang bang control
% There is a bug with the starting point and the sign of the gamma0+tr
gamma_0 = 0;
t_0 = 0;

sat = @(t) min(1,max(0,t));
t_d = @(g_so) abs(2*(g_so-.5));
t_s = @(t_so,g_so) t_so*(1-t_d(g_so));
g_bb = @(t_n,t_so,g_so) sign(g_so-.5)*rn*t_d(g_so)*sat((t_n-t_s(t_so,g_so))/t_d(g_so));


%% Plot bang-bang 2 dof
subplot(2,2,4)
hold on

for nI = 1:nIterations
    % generate 2xn seeds
    o1 = rand(nBranches,1); % t_so
    o2 = rand(nBranches,1); % g_so

    for iO = 1:nBranches
        % calculate trajectory
        tr = g_bb(tn,o1(iO),o2(iO)) + gamma_0;

        if iO == nBranches
            % Plot chosen one
            plot((0:nControlSteps)*deltaT+t_0, tr,'--',"Color",...
                [0 155 119]./255, "LineWidth",1.5)
            plot((nControlSteps:nHorizonSteps)*deltaT+t_0,...
                tr(end)*ones(size(nControlSteps:nHorizonSteps)),'--',"Color",...
                [0 155 119]./255, "LineWidth",1.5)
            plot((0:nUpdateSteps-1)*deltaT+t_0, tr(1:nUpdateSteps),"Color",...
                [108 194 74]./255, "LineWidth",2.5)

            % Add line to mark end of the control horizon.
            plot([nControlSteps, nControlSteps]*deltaT + t_0,...
                [-2.9,2.9]+tr(nControlSteps),...
                "LineWidth",2,"Color",[255,184,28]./255)

            % Assign new starting conditions
            gamma_0 = g_bb(tn(nUpdateSteps),o1(iO),o2(iO)) + gamma_0;
            t_0 = t_0 + (nUpdateSteps-1)*deltaT;
        else
            % Plot disregarded ones
            plot((0:nControlSteps)*deltaT+t_0,tr,"Color", ...
                [194 194 194]./255, "LineWidth",1)
            plot((nControlSteps:nHorizonSteps)*deltaT+t_0,...
                tr(end)*ones(size(nControlSteps:nHorizonSteps)),"Color",...
                [194 194 194]./255, "LineWidth",1)
        end
    end
    % Add line to mark transition between two trajectories.
    plot([t_0, t_0], [-2.9,2.9]+gamma_0,...
        "LineWidth",2,"Color",[237,104,66]./255)
end

hold off
xlabel('time (s)')
ylabel('Orientation (deg)')
xlim([0,450])

%% Part 3 scheduled Bang bang control
gamma_0 = 0;
t_0 = 0;
t_o = 0;

sat = @(t) min(1,max(0,t));
t_d = @(g_so) abs((g_so-.5));
g_bb = @(t_n,g_so) sign(g_so-.5)*rn*t_d(g_so)*sat((t_n-t_o)/t_d(g_so));


%% Bang-bang 1-dof
subplot(2,2,3)
hold on

for nI = 1:nIterations
    % generate 2xn seeds
    o1 = rand(nBranches,1); % g_so

    for iO = 1:nBranches
        % calculate trajectory
        tr = g_bb(tn,o1(iO)) + gamma_0;

        if iO == nBranches
            % Plot chosen one
            plot((0:nControlSteps)*deltaT+t_0, tr,'--',"Color",...
                [108 194 74]./255, "LineWidth",1.5)
            plot((nControlSteps:nHorizonSteps)*deltaT+t_0,...
                tr(end)*ones(size(nControlSteps:nHorizonSteps)),'--',"Color",...
                [108 194 74]./255, "LineWidth",1.5)
            plot((0:nUpdateSteps-1)*deltaT+t_0, tr(1:nUpdateSteps),"Color",...
                [0 155 119]./255, "LineWidth",2.5)

            % Add line to mark end of the control horizon.
            plot([nControlSteps, nControlSteps]*deltaT + t_0,...
                [-2.9,2.9]+tr(nControlSteps),...
                "LineWidth",2,"Color",[255,184,28]./255)

            % Assign new starting conditions
            gamma_0 = g_bb(tn(nUpdateSteps),o1(iO)) + gamma_0;
            t_0 = t_0 + (nUpdateSteps-1)*deltaT;
        else
            % Plot disregarded ones
            plot((0:nControlSteps)*deltaT+t_0,tr,"Color", ...
                [194 194 194]./255, "LineWidth",1)
            plot((nControlSteps:nHorizonSteps)*deltaT+t_0,...
                tr(end)*ones(size(nControlSteps:nHorizonSteps)),"Color",...
                [194 194 194]./255, "LineWidth",1)
        end
    end
    % Add line to mark transition between two trajectories.
    plot([t_0, t_0], [-2.9,2.9]+gamma_0,...
        "LineWidth",2,"Color",[237,104,66]./255)
end

hold off
xlabel('time (s)')
ylabel('Orientation (deg)')
xlim([0,450])


% for iS = 1:4
%     subplot(2,2,iS)
%     ylim([-30,50])
% end
% exportgraphics(gcf,'random_tranjectories.pdf','ContentType','vector')