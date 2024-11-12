[S0,S1,S2] = meshgrid(linspace(-.5,.5),linspace(-1.2,1.2),linspace(-1.2,1.2));

C = true(size(S0));

for x = linspace(0,1)
    C = and(C(:),...
        abs(S0(:)*(1-x)^3*4+S1(:)*(1-x)^2*x*12+S2(:)*(1-x)*x^2*12)<1);
end

%%
figure

[k2,av2] = convhull(S0(C),S1(C),S2(C),'Simplify',true);

trisurf(k2,S0(C),S1(C),S2(C),"edgecolor",'none','facecolor','interp')
axis equal


%scatter3(S0(C),S1(C),S2(C),10,[0,0,0],"filled")
axis equal
xlabel('s_0')
ylabel('s_1')
zlabel('s_2')