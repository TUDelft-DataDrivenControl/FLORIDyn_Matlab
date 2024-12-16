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


function [M,T] = setUpTmpWFAndRun(T,paramFLORIS,Wind)
%SETUPTMPWFANDRUN Summary of this function goes here
%   Detailed explanation goes here

M           = zeros(T.nT,3);
T.weight    = cell(T.nT,1);
T.red_arr   = ones(T.nT,T.nT); % reset
% This should be parallized
for iT = 1:T.nT
    % Interpolate Wind field if needed
    iTWFState = T.States_WF(T.StartI(iT),:);
    if isfield(T,'C_Vel')
        % Weighted OPs are iterated instead of each OP by themselves 
        iTWFState(1) = T.C_Vel(iT,:) * T.States_WF(:,1);
    end
    if isfield(T,'C_Dir')
        % Weighted OPs are iterated instead of each OP by themselves 
        iTWFState(2) = T.C_Dir(iT,:) * T.States_WF(:,2);
    end

    % Single turbine?
    if isempty(T.dep{iT})
        % if not dependent on other turbines, submit as single turbine to
        % only calculate influences from elements like wind shear
        [T_red_arr,~, ~] = runFLORIS(...
            T.posBase(iT,:)+T.posNac(iT,:),... % tmp_Tpos
            iTWFState,...
            T.States_T(T.StartI(iT),:),...
            T.D(iT),...
            paramFLORIS,Wind.Shear);
        M(iT,:) = [T_red_arr, 0, T_red_arr*T.States_WF(T.StartI(iT),1)];
        T.red_arr(iT,iT) = T_red_arr; 
        continue; 
    end
    
    % create wind farm
    tmp_nT   = length(T.dep{iT})+1;
    %   Init states with goal turbine states (last one)
    tmp_Tpos = repmat(T.posBase(iT,:)+T.posNac(iT,:),tmp_nT,1);
    tmp_WF   = repmat(iTWFState,tmp_nT,1);
    tmp_Tst  = repmat(T.States_T(T.StartI(iT),:),tmp_nT,1);
    if T.D(end)>0
        tmp_D = [T.D(T.dep{iT});T.D(iT)];
    else
        % For plotting
        tmp_D = T.D;
    end
    
%     % ================== DEBUG PLOT 1 ========================
%     figure
%     subplot(2,1,1)
%     scatter(tmp_Tpos(end,1),tmp_Tpos(end,2),'filled')
%     hold on
    
    
    % Derive tmp turbine locations and states
    for iiT = 1:(tmp_nT-1)
        % Calculate interpolated OP
        OP1_i = T.intOPs{iT}(iiT,1); % Index OP 1
        OP1_r = T.intOPs{iT}(iiT,2); % Ratio OP 1
        OP2_i = T.intOPs{iT}(iiT,3); % Index OP 2
        OP2_r = T.intOPs{iT}(iiT,4); % Ratio OP 2
        
        % Interpolate location
        OPi_l = ...
            OP1_r * T.States_OP(OP1_i,:) + ...
            OP2_r * T.States_OP(OP2_i,:);
        
        % Apply Interpolated OP temporarily as location
        tmp_Tpos(iiT,:) = OPi_l(1:3);
        
        % Write interpolated turbine states
        tmp_Tst(iiT,:) = OP1_r * T.States_T(OP1_i,:) + ...
            OP2_r * T.States_T(OP2_i,:);
        
        % Write interpolated wind field states
        tmp_WF(iiT,:) = OP1_r * T.States_WF(OP1_i,:) + ...
                OP2_r * T.States_WF(OP2_i,:);
        
        % Potentially overwrite weighted entries
        si = T.StartI(T.dep{iT}(iiT));
        if isfield(T,'C_Vel')
            % Weighted OPs
            tmp_WF(iiT,1) = T.C_Vel(iT,si:si+T.nOP-1)./...
                sum(T.C_Vel(iT,si:si+T.nOP-1)) * ...
                T.States_WF(si:si+T.nOP-1, 1);
        end
        if isfield(T,'C_Dir')
            % Weighted OPs
            tmp_WF(iiT,2) = T.C_Dir(iT,si:si+T.nOP-1)./...
                sum(T.C_Dir(iT,si:si+T.nOP-1)) * ...
                T.States_WF(si:si+T.nOP-1, 2);
        end
        

        % Subtract the wake vector to OP rotated into the real coord.
        % system to find theoretical turbine location
        %   Either use wind dir as OP orientation or original orientation
        if size(tmp_WF,2) == 4
            % Original orientation
            tmp_phi = angSOWFA2world(tmp_WF(iiT,4));
        else
            % Wind dir
            tmp_phi = angSOWFA2world(tmp_WF(iiT,2));
        end
        
        
        tmp_Tpos(iiT,1) = tmp_Tpos(iiT,1) - (...
                cos(tmp_phi) * OPi_l(4) - ... % CHECK IF SIN SIGN IS OK
                sin(tmp_phi) * OPi_l(5));
        
        tmp_Tpos(iiT,2) = tmp_Tpos(iiT,2) - (...
                sin(tmp_phi) * OPi_l(4) + ... % CHECK IF SIN SIGN IS OK
                cos(tmp_phi) * OPi_l(5));
        
        tmp_Tpos(iiT,3) = tmp_Tpos(iiT,3) - OPi_l(6);
        
%         % ================== DEBUG PLOT 2 ========================
%         % OPs of Tii
%         scatter([T.States_OP(OP1_i,1),T.States_OP(OP2_i,1)],...
%             [T.States_OP(OP1_i,2),T.States_OP(OP2_i,2)])
%         % Tii
%         scatter(tmp_Tpos(iiT,1),tmp_Tpos(iiT,2))
%         % Vec from Tii to OP interpolated
%         quiver(tmp_Tpos(iiT,1),tmp_Tpos(iiT,2),...
%             OPi_l(1)-tmp_Tpos(iiT,1),OPi_l(2)-tmp_Tpos(iiT,2),0)
%         % Vec from OP interpolated to turb
%         quiver(OPi_l(1),OPi_l(2),...
%             tmp_Tpos(end,1)-OPi_l(1),tmp_Tpos(end,2)-OPi_l(2),0)
    end
%     % ================== DEBUG PLOT 3 ========================
%     axis equal
%     grid on
%     hold off
%     title('Temp. coordinates')
%     subplot(2,1,2)
%     scatter(T.posBase(iT,1),T.posBase(iT,2),'filled')
%     hold on
%     scatter(T.posBase(T.dep{iT},1),T.posBase(T.dep{iT},2))
%     scatter([T.States_OP(OP1_i,1),T.States_OP(OP2_i,1)],...
%             [T.States_OP(OP1_i,2),T.States_OP(OP2_i,2)])
%     axis equal
%     grid on
%     hold off
%     title('Real coordinates')
    
    %% Run FLORIS sim and get results from last turbine
    [T_red_arr,T_aTI_arr, T_Ueff, T_weight] = ...
        runFLORIS(tmp_Tpos,tmp_WF,tmp_Tst, tmp_D, paramFLORIS,Wind.Shear);
    
    % Wind speed reduction
    T_red        = prod(T_red_arr);
    % Wind speed reduction by individual turbine
    T.red_arr(iT,[T.dep{iT}, iT]) = T_red_arr; 
    % Added turbulence intensity
    T_addedTI    = sum(T_aTI_arr.^2); % <===== CHANGED
    %T_addedTI    = sum(T_aTI_arr); % <===== OLD
    T_addedTI    = sqrt(T_addedTI);

    T.weight{iT} = T_weight;

    if T.D(end)<=0
        % Plotting
        dists   = zeros((tmp_nT-1),1);
        plot_WF = zeros((tmp_nT-1),size(T.States_WF,2));
        plot_OP = zeros((tmp_nT-1),2);
        for iiT = 1:(tmp_nT-1)
            % Calculate interpolated OP
            OP1_i = T.intOPs{iT}(iiT,1); % Index OP 1
            OP1_r = T.intOPs{iT}(iiT,2); % Ratio OP 1
            OP2_i = T.intOPs{iT}(iiT,3); % Index OP 2
            OP2_r = T.intOPs{iT}(iiT,4); % Ratio OP 2
            
            % Interpolate location
            OPi_l = ...
                OP1_r * T.States_OP(OP1_i,:) + ...
                OP2_r * T.States_OP(OP2_i,:);
            plot_OP(iiT,:) = OPi_l(1:2);
            plot_WF(iiT,:) = OP1_r * T.States_WF(OP1_i,:) + ...
                             OP2_r * T.States_WF(OP2_i,:);
            dists(iiT) = sum(sqrt((OPi_l(1:2) - T.posBase(iT,1:2)).^2));
        end
        
        [~,I] = sort(dists);
        
        if length(I) == 1
            Ufree = plot_WF(I(1),1);
            T_Ueff = T_red* Ufree;
        else
            a = plot_OP(I(1),:)';
            b = plot_OP(I(2),:)';
            c = T.posBase(iT,1:2)';
            d = ((b-a)'*(c-a))/((b-a)'*(b-a));
            d = min(max(d,0),1);
            r1 = 1-d;
            r2 = d;
            Ufree = r1*plot_WF(I(1),1) + r2*plot_WF(I(2),1);
            
            
            %         Ufree =  plot_WF(I(1),1)*dists(I(2))/(dists(I(1)) + dists(I(2))) + ...
            %                  plot_WF(I(2),1)*dists(I(1))/(dists(I(1)) + dists(I(2)));
            
            T_Ueff = T_red* Ufree;
        end
    end
    M(iT,:) = [T_red, T_addedTI, T_Ueff];
    
    % Making sure the weights add up to 1 or 0 (if no weights)
    wS = sum(T.weight{iT});
    if wS > 0
        T.weight{iT} = T.weight{iT}./wS;
    else
        T.weight{iT}(:) = 0;
    end
end

end

