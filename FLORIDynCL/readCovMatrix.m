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

function [CovMat, CholSig] = readCovMatrix(CovMatData,nT,csvPrefix)
% READCOVMATRIX interprets the data from the read *Covariance.csv file
%   A single value will be used as the same variance for every turbine, all
%   are independent.
%   If as many values as turbines are given they are used as variance for
%   the turbines which are independent.
%   If nT x nT values are given they are directly used as the covariance
%   matrix
% ======================================================================= %
% INPUT
%   CovMatData  scalar, vector, matrix  Sets the covariance matrix
%   nT          scalar                  Number of turbines
%   csvPrefix   string                  Name of convariance, for debugging
%
% OUTPUT
%   CovMat      matrix                  Identified covariance matrix
%   CholSig     matrix                  cholesky decomposed cov. matrix

if sum(size(CovMatData) == 1) == 2
    CovMat = eye(nT)*CovMatData;
elseif sum(size(CovMatData) == nT)<1
    % wrong dimension
    error([csvPrefix 'Covariance.csv has the wrong size, should '...
        'contain either a 1x' num2str(nT) ' matrix, a ' ...
        num2str(nT) 'x' num2str(nT) ' matrix or a single value.'])
elseif size(CovMatData,1)==1
    CovMat = diag(CovMatData);
else
    CovMat = reshape(CovMatData,nT,nT);
end
CholSig = chol(CovMat);
end