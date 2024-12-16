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

function varargout = EnKF_CalcStateCov(EnStates,nE,varargin)
%ENKF_CALCSTATECOV Calculates the state error covariance matrix based on
%three different methods:
%   Only on the first OP  ->    leads to an not observable model
%   With a time varying C ->    based on the two OPs which influence the 
%                               given turbine
%   Input - Output        ->    Calculates the In-/Output correlation
%
%% Example calls
% Only first OP
% C_xx_Vel = EnKF_CalcStateCov(EnKF.States.Vel,EnKF.nE,T.StartI)
%
% With a time varying C
% C_xx_Vel = ...
%   EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE, 'Method', 'AdaptiveC',...
%   'Interaction', EnKF.Interaction, 'StartI', T.StartI,...
%   'iE', iE, 'nT', T.nT, 'nOP', T.nOP);
% 
% With C based on interaction
% [C_xx_Vel, C] = EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE,...
%   'Method', 'AdaptiveC', 'Interaction', EnKF.Interaction, ...
%   'StartI', T.StartI, 'iE', iE, 'nT', T.nT, 'nOP', T.nOP);
%
% Input/Output correlation
% [C_xy, C_yy] = ...
%   EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE, 'Method', 'InputOutput',...
%   'Output',EnKF.Output.Pow)
%% Defaults
method = 'StartI';
%% Check which method is used
if nargin == 3
    StartI = varargin{1};
elseif nargin > 3
    method = varargin{circshift(strcmp(varargin,'Method'),1)};
    
    if strcmpi(method,'starti')
        StartI = varargin{circshift(strcmp(varargin,'StartI'),1)};
    end
else
    error('EnKF_CalcStateCov: No enough input values')
end

%% Apply calculation
EnStatesMean = mean(EnStates,2);

switch lower(method)
    case 'starti'
        varargout{1} = (EnStatesMean - EnStates) * ...
            (EnStatesMean(StartI) - EnStates(StartI,:))'...
            / (nE - 1);
    case 'adaptivec'
        Interaction = varargin{circshift(strcmp(varargin,'Interaction'),1)};
        iE          = varargin{circshift(strcmp(varargin,'iE'),1)};
        nT          = varargin{circshift(strcmp(varargin,'nT'),1)};
        nOP         = varargin{circshift(strcmp(varargin,'nOP'),1)};
        StartI      = varargin{circshift(strcmp(varargin,'StartI'),1)};
        
        C_xx = zeros(nOP*nT,nT);
        C = zeros(nT,nOP*nT);
        
        for iT = 1:nT
            if sum(Interaction{iE,iT}(:,end))==0
                % No influence
                C(iT,StartI(iT)) = 1;
                C_xx(:,iT) = ...
                    (EnStatesMean - EnStates) * ...
                    (EnStatesMean(StartI(iT)) - ...
                    EnStates(StartI(iT),:))' / (nE - 1);
            else
                % Influenced by other turbines
                for iiT = 1:size(Interaction{iE,iT},1)
                    C(iT,Interaction{iE,iT}(iiT,2)) = ...
                        Interaction{iE,iT}(iiT,3) * ...
                        Interaction{iE,iT}(iiT,6);
                    C(iT,Interaction{iE,iT}(iiT,4)) = ...
                        Interaction{iE,iT}(iiT,5) * ...
                        Interaction{iE,iT}(iiT,6);
                end
                nze = C(iT,:) ~= 0;
                
                C_xx(:,iT) = ...
                    (EnStatesMean - EnStates) * ...
                    ((EnStatesMean(nze) - EnStates(nze,:))' * ...
                    C(iT,nze)') / (nE - 1);
                
                % Slight difference (10^-16) between the two solutions
                % C_xx_Vel(:,iT) = ...
                %     (EnKF.States.Vel - EnKF.States.Vel_mean) * ...
                %     (EnKF.States.Vel - EnKF.States.Vel_mean)'...
                %     / (EnKF.nE - 1) * ...
                %     C(iT,:)';
            end
        end
        varargout{1} = C_xx;
        varargout{2} = C;
    case 'adaptiveweightedc'
        C = varargin{circshift(strcmp(varargin,'C'),1)};
        
        C_xx = (EnStatesMean - EnStates) * (EnStatesMean - EnStates)'/...
            (nE - 1);
        
        varargout{1} = C_xx * C';
    case 'inputoutput'
        OutStates = varargin{circshift(strcmp(varargin,'Output'),1)};
        EnOutputMean = mean(OutStates,2);
        
        % C_xy
        varargout{1} =...
            (EnStatesMean - EnStates) * ...
            (EnOutputMean - OutStates)'/(nE - 1);
        % C_yy
        varargout{2} =...
            (EnOutputMean - OutStates) * ...
            (EnOutputMean - OutStates)'/(nE - 1);
    otherwise
        error(['EnKF_CalcStateCov: State error covariance calculation '...
            'method "' method '" not known.'])
end
end

