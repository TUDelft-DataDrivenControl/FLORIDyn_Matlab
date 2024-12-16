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

function K = EnKF_CalcKalmanGain(method, varargin)
%ENKF_CALCKALMANGAIN calculates the Kalman Gain Matrix based on the
%provided information
% K = C_f,xx * M^T (M C_f,xx M^T + C_ee)^-1
%
% IMPORTANT
%   The code assumes that the state error covariance matrix is provided in
%   the form Cxx = C_f,xx M^T as this reduces the computational cost.
%   Therefore the method has to be chosen fitting to the choice in
%   EnKF_CalcStateCov.m
% 
%%
switch lower(method)
    case 'starti_cxx_cee'
        % Based on State error covariance and output error covariance,
        % based on starting index
        Cxx = varargin{circshift(strcmp(varargin,'Cxx'),1)};
        Cee = varargin{circshift(strcmp(varargin,'Cee'),1)};
        StI = varargin{circshift(strcmp(varargin,'StartI'),1)};
        
        K = Cxx/(Cxx(StI,:) + Cee);
    case 'c_cxx_cee'
        % Based on State error covariance and output error covariance and C
        % matrix
        Cxx = varargin{circshift(strcmp(varargin,'Cxx'),1)};
        Cee = varargin{circshift(strcmp(varargin,'Cee'),1)};
        C   = varargin{circshift(strcmp(varargin,'C'),1)};
        
        K = Cxx/(C * Cxx + Cee);
    case 'c_cxx_cee_loc'
        % Based on State error covariance and output error covariance and C
        % matrix with localization
        Cxx = varargin{circshift(strcmp(varargin,'Cxx'),1)};
        Cee = varargin{circshift(strcmp(varargin,'Cee'),1)};
        C   = varargin{circshift(strcmp(varargin,'C'),1)};
        loc = varargin{circshift(strcmp(varargin,'Loc'),1)};
        
        K = (loc.*Cxx)/(C * (loc.*Cxx) + Cee);
    case 'starti_cxx_cee_loc'
        % Based on State error covariance and output error covariance
        Cxx = varargin{circshift(strcmp(varargin,'Cxx'),1)};
        Cee = varargin{circshift(strcmp(varargin,'Cee'),1)};
        StI = varargin{circshift(strcmp(varargin,'StartI'),1)};
        loc = varargin{circshift(strcmp(varargin,'Loc'),1)};
        
        K = (loc.*Cxx)/(loc(StI,:) * Cxx(StI,:) + Cee);
    case 'cxy_cyy_cee'
        % Based on in-/output correlation
        Cxy = varargin{circshift(strcmp(varargin,'Cxy'),1)};
        Cyy = varargin{circshift(strcmp(varargin,'Cyy'),1)};
        Cee = varargin{circshift(strcmp(varargin,'Cee'),1)};
        
        K = (LocCov.*Cxy)/(Cyy + Cee);
    case 'cxy_cyy_cee_loc'
        % Based on in-/output correlation with localization
        Cxy = varargin{circshift(strcmp(varargin,'Cxy'),1)};
        Cyy = varargin{circshift(strcmp(varargin,'Cyy'),1)};
        Cee = varargin{circshift(strcmp(varargin,'Cee'),1)};
        loc = varargin{circshift(strcmp(varargin,'Loc'),1)};
        
        %K = (loc.*Cxy)/(Cyy + Cee); 
        K = (loc.*Cxy)/(Cyy.* loc(1:200:end,:) + Cee);
end

end

