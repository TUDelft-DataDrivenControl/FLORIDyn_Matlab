function e = error_EnKF(EnKF,TnT,powSOWFA,uSOWFA,e_start)
%ERROR_ENKF Compares the measurements obtained by the EnKF to SOWFA
%reference data - power and effective wind speed. The error quantities aim
%to describe various aspects of the simulation performance.


% Get mean std min max values of turbine measurements
[mM,stdM,minM,maxM] = EnKF_getMeanMeasurements(EnKF);

time = mM.("Time [s]")(1:TnT:end);
nE = 5; % Number of turbine individual error QoIs
e = zeros(6 + TnT*nE,1);
% Calculate errors
%% ======= Power error =======
% Mean squared error turbine power
S_pow = reshape(powSOWFA(:,3),TnT,[])';
S_tim = powSOWFA(1:TnT:end,2);
S_pow_interp    = interp1(S_tim,S_pow,time(time>=e_start));
EnKF_pow        = reshape(mM.("Power generated [MW]"),TnT,[])';
EnKF_pow        = EnKF_pow(time>=e_start,:);
EnKF_pow_std    = reshape(stdM.("Power generated [MW]"),TnT,[])';
EnKF_pow_std    = EnKF_pow_std(time>=e_start,:);

% Mean bias turbine power
e(1) = mean(mean(mean(S_pow_interp,1) - mean(EnKF_pow,1)));
disp(['Out: mean turbine power bias: ' num2str(e(1))])

% Mean absolute turbine power error
e(2) = mean((mean(abs(mean(S_pow_interp,1) - mean(EnKF_pow,1)))));
disp(['Out: mean absolute turbine power error: ' num2str(e(2))])

% Mean log_3 turbine error
% e(2) = mean(mean(log(abs(S_pow_interp - EnKF_pow))./log(3),1));
% disp(['Out: Mean log_3 turbine error: ' num2str(e(2))])

% Mean squared error turbine power
e(3) = mean(mean((S_pow_interp - EnKF_pow).^2,1));
disp(['Out: mean squared error turbine power: ' num2str(e(3))])

% Mean squared error farm power
e(4) = mean((sum(S_pow_interp,2) - sum(EnKF_pow,2)).^2);
disp(['Out: mean squared error farm power: ' num2str(e(4))])

% Bias farm power
e(5) = mean(sum(S_pow_interp,2)) - mean(sum(EnKF_pow,2));
disp(['Out: Bias farm power: ' num2str(e(5))])

% Gaussian weighted mean absolute turbine power error
e(6) = mean(...
    mean(abs((S_pow_interp - EnKF_pow)) .* ...
    (1 - exp(-0.5*(S_pow_interp - EnKF_pow).^2./EnKF_pow_std.^2)),1));
disp(['Out: Gaussian weighted mean absolute turbine power error: ' num2str(e(6))])

% Calculate best correlation and related delay for each turbine, as well as
% other turbine individual values

for iT = 1:TnT
    % correlation and related delay
    [e(7 + nE*(iT-1)), e(8 + nE*(iT-1))] = ...
        corr_error(EnKF_pow(:,iT),S_pow_interp(:,iT),-40:10);
    disp(['Out - T' pad(num2str(iT),2,'left','0') ...
        ': Correlation : ' num2str(e(7 + nE*(iT-1)))])
    disp(['Out - T' pad(num2str(iT),2,'left','0') ...
        ': Delay : ' num2str(e(8 + nE*(iT-1)))])

    % power bias
    e(9 + nE*(iT-1)) = ...
        mean(S_pow_interp(:,iT) - EnKF_pow(:,iT));
    disp(['Out - T' pad(num2str(iT),2,'left','0') ...
        ': Power Bias : ' num2str(e(9 + nE*(iT-1)))])

    % Mean absolute error
    e(10 + nE*(iT-1)) = ...
        mean(abs(S_pow_interp(:,iT) - EnKF_pow(:,iT)));
    disp(['Out - T' pad(num2str(iT),2,'left','0') ...
        ': Mean absolute error : ' num2str(e(10 + nE*(iT-1)))])

    % Mean weighted absolute error
    e(11 + nE*(iT-1)) = ...
        mean( ...
        abs(S_pow_interp(:,iT) - EnKF_pow(:,iT)) .* ...
        (1 - exp(-0.5*(S_pow_interp(:,iT) - EnKF_pow(:,iT)).^2./ ...
        EnKF_pow_std(:,iT).^2))...
        );

    disp(['Out - T' pad(num2str(iT),2,'left','0') ...
        ': Mean weighted absolute error : ' num2str(e(11 + nE*(iT-1)))])
end
%% ======= Wind speed error =======
% Mean squared error turbine wind speed

% Bias turbine wind speed

end

function [c,d] = corr_error(enkf,ref,d_test)
% Calculate the correlation between the output of the enkf and the
% reference data. Goal is to determine the adequate delay bewteen the two
% and to strengthen the correlation
%
% INPUTS
%   enkf    Values estimated by the ensemble kalman filter
%   ref     Reference values
%   d_test  Test time step differences
%               e.g. for +- 5 time steps -> d_test = -5:5
%
% OUTPUTS
%   c       Maximum correlation value
%   d       Respective delay
%               d < 0 : enkf data is too early / too fast
%               d > 0 : enkf data is too late / too slow
%%
corr_value = zeros(size(d_test));

if max(abs(d_test))>length(enkf)
    disp(['ERROR Correlation calculation: Number of delays:' ...
        num2str(length(d_test)) ', number of values: ' ...
        num2str(length(enkf)) '.'])

    if length(enkf)>1
        % Too long range was requested, reducing it.
        disp('Restructuring delays to be tested ...')
        d_test = (1:length(enkf)-1) - ceil(length(enkf)/2);
    end
end

try
    for iD = 1:length(d_test)
        if abs(d_test(iD))>length(enkf)+2; corr_value(iD) = -2; continue; end
        if d_test(iD)<0
            % Start EnKF data later and end reference data earlier
            corr_value(iD) = ...
                corr(enkf(-d_test(iD)+1:end), ref(1:end+d_test(iD)));
        else
            % End EnKF data earlier and start reference data later
            corr_value(iD) = ...
                corr(enkf(1:end-d_test(iD)), ref(d_test(iD)+1:end));
        end
    end
    % Set outputs
    [c,I]   = max(corr_value);
    d       = d_test(I);
catch
    % Error during calculation
    disp('ERROR Correlation calculation: resorted to c = 0')
    c = 0;
    d = 0;
end

%

%% Debug
debug = false;
if debug
    figure
    subplot(3,1,1)
    plot(d_test,corr_value,'k','LineWidth',2)
    hold on
    scatter(d,c,20,"cyan")
    
    subplot(3,1,2)
    plot(enkf,'--','LineWidth',1.5,'Color',[.7,.7,.7])
    hold on
    plot(ref,'-','LineWidth',2,'Color',"red")
    if d<0
        plot(-d_test(iD):length(enkf), ...
            enkf(-d_test(iD)+1:end),'-k','LineWidth',2)
    else
        plot(1:length(enkf)-d_test(iD), ...
            enkf(1:end-d_test(iD)),'-k','LineWidth',2)
    end
    hold off
    legend('EnKF', 'Ref','Corrected')
end

end