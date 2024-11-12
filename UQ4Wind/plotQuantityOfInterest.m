function plotQuantityOfInterest(time,data,TnT,figNr,dataName,prog)
%PLOTQUANTITYOFINTEREST Summary of this function goes here
%   Detailed explanation goes here

f = figure(figNr);
c = TUDelft_Seq_wyr(1000);
col = c(min(round(prog*1000),1000),:);
lay = getLayout(TnT);
for iT = 1:TnT
    subplot(lay(1),lay(2),iT)
    hold on
    plot(time(iT:TnT:end),data(iT:TnT:end),'Color',col)
    hold off
    grid on
    xlabel('Time [s]')
    ylabel(dataName)
    xlim([time(1),time(end)]);
end
end

