function r = EnKF_GetAveForeignReduction(M,nT)

r = zeros(nT,1);
for iE = 1:length(M)
    r = r + M{iE}.("Foreign Reduction [%]")(end-nT+1:end);
end
r = r./length(M);

end