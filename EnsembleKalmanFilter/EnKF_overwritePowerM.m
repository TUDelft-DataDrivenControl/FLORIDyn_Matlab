function EnKF = EnKF_overwritePowerM(EnKF, TnT, P)

for iE = 1:EnKF.nE
    EnKF.M{iE}.("Power generated [MW]")(end-TnT+1:end) = P(:,iE);
end

end