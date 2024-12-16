costFun = readmatrix('WindDir.csv');
costFun(:,2) = costFun(:,2) - costFun(1,2);
costFun(:,2) = costFun(:,2)./costFun(end,2);

Ustart  = readmatrix('probe_T0_start_U.csv');
U = zeros(size(Ustart,1),10);
U(:,1) = Ustart(:,1);

for iT = 0:8
    Ustart    = readmatrix(['probe_T' num2str(iT) '_start_U.csv']);
    Uend      = readmatrix(['probe_T' num2str(iT) '_end_U.csv']);
    Uint      = interp1(costFun(:,1),costFun(:,2),Ustart(:,1));
    U(:,iT+2) = Ustart(:,2).*(1-Uint) + Uend(:,2).*Uint;
end

writematrix(U,'U.csv');

%%
figure
hold on
for i = 1:9
    plot(U(:,1)-U(1,1),U(:,1+i));
end
hold off
grid on
