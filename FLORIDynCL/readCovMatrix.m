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