figure


subplot(1,3,1)
contourf(C_xx_Dir,'EdgeColor','none');
axis equal
colormap(viridis(1000));
title('Original Convariance')
set(gca,'YDir','reverse');

subplot(1,3,2)
contourf(GaspariAndCohn1999(EnKF.Dir.cutOffLength,distOPs),'EdgeColor','none');
axis equal
colormap(viridis(1000));
title('Localization Convariance')
set(gca,'YDir','reverse');

subplot(1,3,3)
contourf(C_xx_Dir,'EdgeColor','none');
axis equal
colormap(viridis(1000));
title('Combined Convariance')
set(gca,'YDir','reverse');


%%
figure

for i=1:9
    subplot(3,3,i)
    hold on
    for iE=1:EnKF.nE
        plot(EnKF.States.Vel(T.StartI(i):T.StartI(i)+T.nOP-1,iE))
    end
    grid on
end