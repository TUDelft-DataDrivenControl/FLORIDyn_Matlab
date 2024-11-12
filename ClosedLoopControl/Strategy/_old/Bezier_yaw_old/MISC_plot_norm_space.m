rate_lim = 1;

figure
hold on
for ix = 0:0.1:1
    y = BC_3_2dof(rate_lim,[ones(1,100)*ix;linspace(0,1,100)]');
    
    if or(ix == 0, ix == 1)
        plot(y(:,1),y(:,2),'Color',[12,35,64]/255,'LineWidth',2);
    else
        plot(y(:,1),y(:,2),'Color',[237,104,66]/255,'LineWidth',2);
    end
end

for iy = 0:0.1:1
    y = BC_3_2dof(rate_lim,[linspace(0,1,100);ones(1,100)*iy]');
    
    if or(iy == 0, iy == 1)
        plot(y(:,1),y(:,2),'Color',[12,35,64]/255,'LineWidth',2);
    else
        plot(y(:,1),y(:,2),'Color',[237,104,66]/255,'LineWidth',2);
    end
end

% plot([-1,1]*rate_lim/3,[-2,2]*rate_lim/3,'--',...
%     'Color',[165,0,52]/255,'LineWidth',2)
% 
% plot([-1,1]*rate_lim/3,[-1,1]*rate_lim/3,'--',...
%     'Color',[165,0,52]/255,'LineWidth',2)

axis equal
ylabel('\gamma_2')
xlabel('\gamma_1')
grid on
xlim([-0.4,0.4])
ylim([-0.8,0.8])

% exportgraphics(gcf,'name.pdf','ContentType','vector')