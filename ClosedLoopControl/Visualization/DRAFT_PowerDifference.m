figure(1)
subplot(2,1,1)
hold on

timePow = Pwr_out(41:290,1);
refPow  = Pwr_out(41:290,end);
genPow  = sum(Pwr_out(41:290,2:end-1),2);

plot(timePow(refPow<genPow),(genPow(refPow<genPow) - refPow(refPow<genPow))./refPow(refPow<genPow),...
    '+','Color',[165, 0, 52]./255,'LineWidth',2)
plot(timePow(refPow>genPow),(genPow(refPow>genPow) - refPow(refPow>genPow))./refPow(refPow>genPow),...
    'o','Color',[165, 0, 52]./255,'LineWidth',2)

hold off



subplot(2,1,2)
hold on

csum_Gen = cumsum(genPow)./cumsum(refPow);

plot(timePow,csum_Gen,'Color',[165, 0, 52]./255,'LineWidth',2)
%plot(timePow(csum>0),csum(csum>0),'Color',[165, 0, 52]./255,'LineWidth',2)
%plot(timePow(csum<0),csum(csum<0),'Color',[165, 0, 52]./255,'LineWidth',2)

hold off