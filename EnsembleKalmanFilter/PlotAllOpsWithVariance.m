stdOPx = std(EnKF.States_OP(:,1:6:end),[],2);
aveOPx = mean(EnKF.States_OP(:,1:6:end),2);

stdOPy = std(EnKF.States_OP(:,2:6:end),[],2);
aveOPy = mean(EnKF.States_OP(:,2:6:end),2);

Colors = [...
    12,35,64 
    255,184,28 
    108,194,74
    0,118,194
    237,104,66 
    0,155,119 
    224,60,49
    111,29,119 
    0,184,200 
    239,96,163 
	165,0,52]./255;
figure
hold on
for i=1:length(stdOPx)
    plot([aveOPx(i)+2*stdOPx(i),aveOPx(i)-2*stdOPx(i)],...
        [aveOPy(i),aveOPy(i)],'-|',...
        'Color',Colors(mod(i,11)+1,:))
    
    plot([aveOPx(i),aveOPx(i)],...
        [aveOPy(i)+2*stdOPy(i),aveOPy(i)-2*stdOPy(i)],'-_',...
        'Color',Colors(mod(i,11)+1,:))
    
end
axis equal
grid on
hold off