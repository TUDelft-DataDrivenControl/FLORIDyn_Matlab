function time_matrix = calcTimes(T, Uinfty)
%CALCTIMES calculates the time between the turbines based on the free wind
%speed and direction

% http://dx.doi.org/10.1016/j.renene.2017.06.065 -> 0.7
% parameter tuning -> 0.7396
U_adv = 0.7396*Uinfty;

% Time matrix lists the time it takes for the actuation signal at T_col to
% reach T_row
time_matrix = zeros(T.nT, T.nT);

for iT = 1:T.nT
    dist    = T.posBase(:,1:2)' - T.posBase(iT,1:2)';
    phi     = angSOWFA2world(T.States_WF(T.StartI(iT),2));
    dist    = [cos(phi), sin(phi); -sin(phi), cos(phi)] * dist;

    time_matrix(iT,:) = dist(1,:)./U_adv;
    
    % Remove entries that are upstream (-> -1), too far away cross stream
    % (-> -2)
    time_matrix(iT,time_matrix(iT,:)<0) = -1;
    time_matrix(iT,and(time_matrix(iT,:)>0, abs(dist(2,:)./T.D(iT))>2)) = -2;
end


end

