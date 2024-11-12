y = BC_3_2dof(rate_lim,[.2,.5; .2,1; .1,.3;.9,.8]);

[t,yaw] = BC_3_2dof_trajectory([1,1,1,1]', 0.01, y, [0;0;0;0]);

figure

plot(t,yaw(:,1),t,yaw(:,2),t,yaw(:,3),t,yaw(:,4),'LineWidth',2)
legend('0.2,0.5','0.2,1','0.1,0.3','0.9,0.8')
grid on
xlabel('t_n')