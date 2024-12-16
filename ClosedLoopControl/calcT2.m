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

function [T_star,combined,ah,ph,solve_order,t2o, t2g, opt_ph] = calcT2(T, Uinfty, ah, dof, deltaT)
%CALCT2 This function derives a number of optimisation problems from the
%current free wind speed and turbine position. This version of the function
%also derives the necesssary prediction horizon.
%Different to calcT is in this version that dof and ah are decoupled. This
%is necessary if e.g. 2 dof describe the control action across the entire
%ah.
%
% INPUTS
% T       : (struct)    Turbine data and OP states
% Uinfty  : (float)     Free wind speed [m/s]
% ah      : (int)       Action horizon in time steps [-]
% deltaT  : (float)     discretization time step
%
% ----------------------------------------------------------------------- %
% OUTPUTS
% T_star  : (nO x nT*ph matrix) maps the power generated from all turbines
%           over the prediction horizon to nO optimisation const functions
% combined: (nO x nT*ah matrix) shows which action horizon steps are
%           optimised in which cost function
% ah      : (int) Action horizon in time steps [-]
% ph      : (int) Prediction horizon in time steps [-]
% solve_order : (1xnC cell array) collection of optimisation problems that
%               can be solved independently and if they depend on
%               oneanother then they are ordered by solve order
% t2o     : (nO x nT matrix) which turbines are part of which optimization 
% t2g     : (nG x nT matrix) depicting which turbines need to simulated
%               together
% opt_ph  : (nO x 1 vector) with the minimal prediction horizon for each
%               optimization problem.
% ----------------------------------------------------------------------- %
% nC: Number of clusters (disconnected cost functions & optimisation)
% nT: Number of turbines
% nO: Number of Optimisation problems
% ----------------------------------------------------------------------- %

%% Retrieve time delays between turbines 
% -1 if upstream, -2 if too far away crossstream
time_matrix = calcTimes(T, Uinfty);
time_matrix = floor(time_matrix./deltaT);
%time_matrix(time_matrix>0) = time_matrix(time_matrix>0)-1;
ph = max(time_matrix,[],"all") + ah + 0;

debug = true;

%% Create T_star
T_star = zeros(T.nT, T.nT*ph);
% Fill T_star
for iT = 1:T.nT
    % Influence of the turbine onto itself
    T_star(iT,(1:ah)+(iT-1)*ph) = 1;
    
    if sum(time_matrix(iT,:)>0) == 0;continue;end   % No downstream turbine

    % Affected downstream turbines with delay
    ds_T = find(time_matrix(iT,:)>0);
    for iiT = 1:length(ds_T)
        offset = time_matrix(iT,ds_T(iiT));
        T_star(iT,(ds_T(iiT)-1)*ph+offset+(1:ah)) = 1;
    end
end



if debug
    figure
    subplot(3,1,1)
    imagesc(T_star)
    title('T^*_1')
end
%% Transformation 1
% Combine rows with at least one common 1 entry
combining   = true;
current_row = 1;
check_row   = 2;
combined    = eye(size(T_star,1));

while combining
    if sum(and(T_star(current_row,:),T_star(check_row,:)))>0
        % Have at least one common entry -> combine
        T_star(current_row, :) = ...
            or(T_star(current_row, :),T_star(check_row, :));
        combined(current_row, :) = ...
            or(combined(current_row, :),combined(check_row, :));
        T_star(check_row,:)     = [];
        combined(check_row,:)   = [];

        % Need to go back to see if previously ignored rows now need to be
        % considered
        check_row = current_row + 1;
    else
        % Have no common rows -> move on
        check_row = check_row + 1;
    end
    
    % Check if check_row is outside of array
    if check_row>size(T_star,1)
        current_row = current_row + 1;
        check_row = current_row + 1;
    end

    if current_row >= size(T_star,1)
        combining = false;
    end
end

if debug
    subplot(3,1,2)
    imagesc(T_star)
    title('T^*_2')
end
%% Transformation 2
% Combine rows that follow one another
combining   = true;
current_row = 1;
check_row   = 2;

% Create shifted row mask
mask_one = repmat((1:ph-1),1,T.nT);
mask_two = repmat((2:ph),1,T.nT);

add_off = repmat((0:T.nT-1), ph-1, 1);
mask_one = mask_one + add_off(:)'*ph;
mask_two = mask_two + add_off(:)'*ph;

compatible = eye(size(T_star,1));

% Find out which rows are compatible
while combining
    combine = true;
    if check_row>size(T_star,1) && current_row >= size(T_star,1); break; end
    for iO = 1:size(T_star,2)-1
        if (T_star(current_row,iO)==1 && T_star(current_row,iO+1)==0 && T_star(check_row,iO+1)==0)
            % Check if at an end of a true row and if it continues in the
            % checked row - if not don't combine
            combine = false;
        end

        if combine == false
            break;
        end
    end

    if combine
        compatible(current_row, check_row) = true;
    end

    % if isequal(T_star(current_row,mask_one),...
    %         T_star(check_row,mask_two))
    %     % Check row does follow current row -> combine
    %     compatible(current_row, check_row) = 1;
    % end

    check_row = check_row + 1;

    % Check if check_row is outside of array
    if check_row>size(T_star,1)
        current_row = current_row + 1;
        check_row = current_row + 1;
    end

    if current_row >= size(T_star,1)
        combining = false;
    end
end

combining   = true;
current_row = 1;
check_row   = 2;

while combining
    if check_row>size(T_star,1) && current_row >= size(T_star,1); break; end

    if sum(and(compatible(current_row,:),compatible(check_row,:)))>0
        % Have at least one common entry -> combine
        compatible(current_row, :) = ...
            or(compatible(current_row, :),compatible(check_row, :));
        T_star(current_row, :) = ...
            or(T_star(current_row, :),T_star(check_row, :));
        combined(current_row, :) = ...
            or(combined(current_row, :),combined(check_row, :));

        T_star(check_row,:)      = [];
        combined(check_row,:)    = [];
        compatible(check_row,:) = [];
    else
        % Have no common rows -> move on
        check_row = check_row + 1;
    end
    
    % Check if check_row is outside of array
    if check_row>size(T_star,1)
        current_row = current_row + 1;
        check_row = current_row + 1;
    end

    if current_row >= size(T_star,1)
        combining = false;
    end
end

if debug
    subplot(3,1,3)
    imagesc(T_star)
    title('T')
end

%% Determine dependencies
opt_dep = eye(size(T_star,1));

% Go through the blocks of each turbine and determine which opt problem
% follows which other one
for iT = 1:T.nT
    block = T_star(:,(iT-1)*ph+1:iT*ph);
    iO = 0;
    for it = 1:ph
        newO = find(block(:,it)==1);
        if ~isempty(newO)
            if iO ~= 0
                opt_dep(newO,iO) = 1;
            end
            iO = newO;
        end
    end
end

if debug
    figure; 
    imagesc(opt_dep)
    title("Optimisation problem dependency")
end

%% Determine solve order
solve_order = cell(0);
ordered = [];
iO = 1;
opt_dep = opt_dep - eye(size(T_star,1));
visited = false(size(T_star,1),1);

for iOpt = 1:size(opt_dep,1)
    if sum(ordered==iOpt)>0; continue; end
    if sum(opt_dep(:,iOpt)>0)>0; continue; end
    [solve_order{iO}, visited] = recursiveSolve(opt_dep, iOpt, visited);
    ordered = [ordered, solve_order{iO}];
    iO = iO + 1;
end

%% Calculate which turbine appears in which optimization
t2o         = zeros(size(T_star,1),T.nT);
t_ph        = zeros(T.nT*ph,1);
t_ph(1:ph)  = 1;
for iT = 1:T.nT
    tmp_t       = T_star*t_ph;
    t2o(:,iT)   = tmp_t>0;
    t_ph        = circshift(t_ph,ph);
end
clear t_ph tmp_t

%% Convey groups of turbines that need to be simulated together
t2g = zeros(length(solve_order),T.nT);
for iG = 1:length(solve_order)
    t2g(iG,:) = sum(t2o(solve_order{iG},:),1)>0;
end

%% Derive prediction horizons for each problem
%time_matrix
opt_ph = zeros(size(t2o,1),1);

for iO = 1:size(t2o,1)
    opt_ph(iO) = max(time_matrix(t2o(iO,:)>0,t2o(iO,:)>0),[],"all") + ah;
end

end

function [order, visited] = recursiveSolve(opt_dep, iO, visited)
% Recursive function to determine the solve order of the optimisation 
% problems
if visited(iO)
    warning(['Optimization problem ' num2str(iO) ' was already visited!!'])
    writematrix(opt_dep,['ERROR_Opt' num2str(iO) 'was_already_visited.csv'])
    order = [];
    return
end

visited(iO) = true;

if ~any(opt_dep(iO,:))
    % No dependency, start of order
    order = iO;
else
    % Existing dependency
    order = [];
    ind_opt = find(opt_dep(iO,:)==1);
    for i = 1:length(ind_opt) % If length(ind_opt)>1 there are brances
        [tmp_order, visited] = recursiveSolve(opt_dep, ind_opt(i), visited);
        order = [tmp_order, order];
    end
    % Add itself to the order
    order = [order, iO];
end

end

