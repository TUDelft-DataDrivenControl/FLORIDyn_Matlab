[X,Y,Z] = meshgrid(linspace(0,1,size(data,1)),linspace(0,1,size(data,2)),9:54);

data_n = data;
for iS = 1:size(data,3)
    max_v = max(data_n(:,:,iS),[],'all');
    min_v = min(data_n(:,:,iS),[],'all');

    if max_v == min_v; continue; end
    data_n(:,:,iS) = (data_n(:,:,iS)-min_v)./(max_v-min_v);
end

% Transparency Map
amap = zeros(20,1);
amap(1:5)=1;
amap(6:10)=linspace(1,0.05,5);
amap(11:end)=.05;

figure
h = slice(X,Y,Z,data_n(:,:,9:end),[],[],[1:1/3:54]);
%h = slice(X,Y,Z,data_n(:,:,9:end),[0:.05:1],[0:.05:1],[]);
set(h,'EdgeColor','none',...
'FaceColor','interp')
alpha('color')

%alphamap('rampdown')
alphamap(amap)

colormap(viridis(10))
c = colorbar;
c.Label.String = 'J_{norm} (-)';
zlabel('Time (25s steps)')
xlabel('Optimisation param. T0')
ylabel('Optimisation param. T1')