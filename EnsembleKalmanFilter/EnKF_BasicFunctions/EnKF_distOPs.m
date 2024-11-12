function distOPs = EnKF_distOPs(OPs_tmp,StartI)
%ENKF_DISTOPS calculates the distance between the OPs. Calculation spread
%out to save computational cost
%% ============================================================= %%
distOPs = (OPs_tmp(:,1) - OPs_tmp(StartI,1)').^2;
distOPs = distOPs + (OPs_tmp(:,2) - OPs_tmp(StartI,2)').^2;
distOPs = distOPs + (OPs_tmp(:,3) - OPs_tmp(StartI,3)').^2;
distOPs = sqrt(distOPs);
end

