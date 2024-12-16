% Copyright (C) <2024>, M Becker
%
% List of the contributors to the development of FLORIDyn: see LICENSE file.
% Description and complete License: see LICENSE file.
	
% This program (FLORIDyn) is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program (see COPYING file).  If not, see <https://www.gnu.org/licenses/>.
% ======================================================================= %
% Updated: 16. Dez. 2024, M. Becker
% ======================================================================= %

function [mM,stdM,minM,maxM] = EnKF_getMeanMeasurements(EnKF)
%ENKF_GETMEANMEASUREMENTS Returns combined measurements of the ensembles as
%mean, standard deviation, min & max

Time    = EnKF.M{1}.('Time [s]');
M_comb  = cat(1,EnKF.M{:});

mM   = zeros(size(EnKF.M{1}));
stdM = zeros(size(EnKF.M{1}));
minM = zeros(size(EnKF.M{1}));
maxM = zeros(size(EnKF.M{1}));
mM(:,1)     = Time;
stdM(:,1)   = Time;
minM(:,1)   = Time;
maxM(:,1)   = Time;


Data_raw = reshape(...
    M_comb.("Foreign Reduction [%]"),[],EnKF.nE);
mM(:,2)     = mean(Data_raw,2);
stdM(:,2)   = std(Data_raw,[],2);
minM(:,2)   = min(Data_raw,[],2);
maxM(:,2)   = max(Data_raw,[],2);

Data_raw = reshape(...
    M_comb.("Added turbulence [%]"),[],EnKF.nE);
mM(:,3)     = mean(Data_raw,2);
stdM(:,3)   = std(Data_raw,[],2);
minM(:,3)   = min(Data_raw,[],2);
maxM(:,3)   = max(Data_raw,[],2);

Data_raw = reshape(...
    M_comb.("Effective Wind Speed [ms^-1]"),[],EnKF.nE);
mM(:,4)     = mean(Data_raw,2);
stdM(:,4)   = std(Data_raw,[],2);
minM(:,4)   = min(Data_raw,[],2);
maxM(:,4)   = max(Data_raw,[],2);

Data_raw = reshape(...
    M_comb.("Free Wind Speed [ms^-1]"),[],EnKF.nE);
mM(:,5)     = mean(Data_raw,2);
stdM(:,5)   = std(Data_raw,[],2);
minM(:,5)   = min(Data_raw,[],2);
maxM(:,5)   = max(Data_raw,[],2);

Data_raw = reshape(...
    M_comb.("Power generated [MW]"),[],EnKF.nE);
mM(:,6)     = mean(Data_raw,2);
stdM(:,6)   = std(Data_raw,[],2);
minM(:,6)   = min(Data_raw,[],2);
maxM(:,6)   = max(Data_raw,[],2);

mM = array2table(mM,...
    'VariableNames',{'Time [s]','Foreign Reduction [%]',...
    'Added turbulence [%]','Effective Wind Speed [ms^-1]',...
    'Free Wind Speed [ms^-1]','Power generated [MW]'});
stdM = array2table(stdM,...
    'VariableNames',{'Time [s]','Foreign Reduction [%]',...
    'Added turbulence [%]','Effective Wind Speed [ms^-1]',...
    'Free Wind Speed [ms^-1]','Power generated [MW]'});
minM = array2table(minM,...
    'VariableNames',{'Time [s]','Foreign Reduction [%]',...
    'Added turbulence [%]','Effective Wind Speed [ms^-1]',...
    'Free Wind Speed [ms^-1]','Power generated [MW]'});
maxM = array2table(maxM,...
    'VariableNames',{'Time [s]','Foreign Reduction [%]',...
    'Added turbulence [%]','Effective Wind Speed [ms^-1]',...
    'Free Wind Speed [ms^-1]','Power generated [MW]'});

end

