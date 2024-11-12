figure
plot(mean(m2(1:9:end,:)))

hold on
plot(m1(1,:))

%%
figure
for iE = 1:EnKF.nE
    for iT = 1:9
        VisTime = EnKF.M{iE}.("Time [s]");
        VisPow  = EnKF.M{iE}.("Foreign Reduction [%]");
        subplot(3,3,iT)
        hold on
        plot(VisTime(iT:T.nT:end),VisPow(iT:T.nT:end),...
            'Color',EnKF.Colors(mod(iE,11)+1,:))
        hold off
        grid on
    end
end