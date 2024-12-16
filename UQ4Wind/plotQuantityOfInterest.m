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

