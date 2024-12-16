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

function U = GetFlowField(T,Wind,Vis,paramFLORIS)
%PLOTFLOWFIELD Generate full flow field plot
% exportgraphics(gcf,'9T_flow.pdf','ContentType','vector')
%% Preallocate field
nM = countMeasurements(Vis.FlowField.Data);
nM = 3;
fieldLims = Vis.FlowField.Lims;
fieldRes  = Vis.FlowField.Res;
xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);

[X,Y] = meshgrid(xAx,yAx);

clear fieldRes 

%% Get data

if Vis.FlowField.Plot.parallel
    Z = getMeasurementsPar2(X,Y,nM,T.posNac(1,3),T,paramFLORIS,Wind,Vis);
else
    Z = getMeasurements(X,Y,nM,T.posNac(1,3),T,paramFLORIS,Wind,Vis);
end

U = Z(:,:,3);

end
%exportgraphics(gcf,'9T_HorFlow_t300.pdf','ContentType','vector')
%%
function Z = getMeasurements(X,Y,nM,zh,T,paramFLORIS,Wind,Vis)
% Wrapper function to disguise the grid points as turbines with one rotor
% point and experience almost the same calculations as the rotor points in
% the simulation.
%   Single thread version
sizeX = size(X);
Z = zeros(sizeX(1),sizeX(2),nM);
%parfor iGP = 1:length(X(:))
for iGP = 1:length(X(:))
    xGP = X(iGP);
    yGP = Y(iGP);
    
    GPdep = {1:T.nT};


    % Create T equivalent to get interpolated OPs
    GP.dep          = cell(T.nT,1);
    GP.dep(1)       = GPdep(1);
    GP.StartI       = T.StartI;
    GP.nT           = 1;
    GP.States_OP    = T.States_OP;
    GP.posBase      = [xGP,yGP,0];
    GP.nOP          = T.nOP;
    GP.intOPs = interpolateOPs(GP);
    
    GP.posNac = [0,0,zh];
    GP.States_WF = T.States_WF;
    GP.States_T = T.States_T;
    GP.D = [T.D;0];
    
    [tmpM,~] = setUpTmpWFAndRun(GP,paramFLORIS,Wind);
    
    [rw,cl] = ind2sub(size(X),iGP);
    Z(rw,cl,1:3) = tmpM;
end
end

function Z = getMeasurementsPar2(X,Y,nM,zh,T,paramFLORIS,Wind,Vis)
% Wrapper function to disguise the grid points as turbines with one rotor
% point and experience almost the same calculations as the rotor points in
% the simulation.
%   Multi thread version, requires parallel computing toolbox
%   For one plot, the start-up time of the parallel toolbox can be longer
%   than the saved computational time.
Zred  = zeros(size(X));
Zeff  = zeros(size(X));
ZambI = zeros(size(X));

%parfor iGP = 1:length(X(:))
parfor iGP = 1:length(X(:))
    xGP = X(iGP);
    yGP = Y(iGP);
    
    GPdep = {1:T.nT};

    GP = struct();
    % Create T equivalent to get interpolated OPs
    GP.dep          = cell(T.nT,1);
    GP.dep(1)       = GPdep(1);
    GP.StartI       = T.StartI;
    GP.nT           = 1;
    GP.States_OP    = T.States_OP;
    GP.posBase      = [xGP,yGP,0];
    GP.nOP          = T.nOP;
    GP.intOPs = interpolateOPs(GP);
    
    GP.posNac = [0,0,zh];
    GP.States_WF = T.States_WF;
    GP.States_T = T.States_T;
    GP.D = [T.D;0];
    
    [tmpM,~] = setUpTmpWFAndRun(GP,paramFLORIS,Wind);
    
    Zred(iGP)  = tmpM(1);
    ZambI(iGP) = tmpM(2);
    Zeff(iGP)  = tmpM(3);
end

Z = cat(3,Zred,ZambI,Zeff);
end

function Z = getMeasurementsDelaunay(T)
DT = delaunay(T.States_OP(:,1),T.States_OP(:,2));
trisurf(DT,T.States_OP(:,1),T.States_OP(:,2),T.States_WF(:,1),...
    'EdgeColor','none','facecolor', 'interp')
end
%%
function Z = Wrap_Z(Z,tmpM,sizeX)
[rw,cl] = ind2sub(sizeX,iGP);
    tmpM1 = sub2ind(sizeZ,rw,cl,1);
    tmpM2 = sub2ind(sizeZ,rw,cl,2);
    tmpM3 = sub2ind(sizeZ,rw,cl,3);
    
    Z(tmpM1) = tmpM(1);
end
%%
function n = countMeasurements(VisFlowFieldData)
n = 0;
if VisFlowFieldData.WindSpeedFree
    n = n+1;
end
if VisFlowFieldData.WindSpeedEff
    n = n+1;
end
if VisFlowFieldData.Reduction
    n = n+1;
end
if VisFlowFieldData.AddedTurb
    n = n+1;
end
if VisFlowFieldData.WindDirection
    n = n+1;
end
end
%exportgraphics(gcf,'name.pdf','ContentType','vector')
%c = colorbar;
%c.Limits = [0,14];
%set(gca, 'Clim', [0, 14])
%
%