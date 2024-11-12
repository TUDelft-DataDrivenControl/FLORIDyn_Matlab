figure

for iE = 1:EnKF.nE
    for iT = 1:9
        % Visualization
        % TODO
        VisTime = EnKF.M{iE}.("Time [s]");
        VisPow  = EnKF.M{iE}.("Power generated [MW]");
        subplot(3,3,iT)
        hold on
        plot(VisTime(iT:T.nT:end),VisPow(iT:T.nT:end),...
            'Color',EnKF.Colors(mod(iE,11)+1,:))
        hold off
        grid on
    end
end