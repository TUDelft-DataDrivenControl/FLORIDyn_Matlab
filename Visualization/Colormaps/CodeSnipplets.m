%% Before plotting
figure
if dark_mode; set(gcf,'Color','k'); end

%% Isosurface plotting example
s = isosurface(X,Y,Z,U,8);
p = patch(s);
set(p,'FaceColor',[1,1,1]);
set(p,'FaceAlpha',0.1);

%% After plotting
ax = gca;
ax.YColor = 'w';
ax.XColor = 'w';
ax.Color = 'k';
ax.GridColor = 'w';