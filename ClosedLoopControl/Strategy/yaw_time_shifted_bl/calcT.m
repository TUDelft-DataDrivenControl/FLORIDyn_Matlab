function [T_star,combined,ah,ph,solve_order] = calcT(T, Uinfty, ah, deltaT)
%CALCT This function derives a number of optimisation problems from the
%current free wind speed and turbine position. This version of the function
%also derives the necesssary prediction horizon
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
%
% ----------------------------------------------------------------------- %
% nC: Number of clusters (disconnected cost functions & optimisation)
% nT: Number of turbines
% nO: Number of Optimisation problems
% ----------------------------------------------------------------------- %

%% Retrieve time delays between turbines 
% -1 if upstream, -2 if too far away crossstream
time_matrix = calcTimes(T, Uinfty);
time_matrix = floor(time_matrix./deltaT);
ph = max(time_matrix,[],"all") + ah +0;

debug = false;

%% Create T_star
T_star = zeros(T.nT*ah, T.nT*ph);
% Fill T_star
for iT = 1:T.nT
    T_star((1:ah)+(iT-1)*ah,(1:ah)+(iT-1)*ph) = eye(ah);
    if sum(time_matrix(iT,:)>0) == 0;continue;end   % No downstream turbine
    ds_T = find(time_matrix(iT,:)>0);
    for iiT = 1:length(ds_T)
        offset = time_matrix(iT,ds_T(iiT));
        for iAH = 1:ah
            T_star((iT-1)*ah + iAH,(ds_T(iiT)-1)*ph+offset+iAH) = 1;
        end
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
        T_star(check_row,:) = [];
        combined(check_row,:) = [];
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
    if isequal(T_star(current_row,mask_one),...
            T_star(check_row,mask_two))
        % Check row does follow current row -> combine
        compatible(current_row, check_row) = 1;
    end

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

for iOpt = 1:size(opt_dep,1)
    if sum(ordered==iOpt)>0; continue; end
    if sum(opt_dep(:,iOpt)>0)>0; continue; end
    solve_order{iO} = recursiveSolve(opt_dep, iOpt);
    ordered = [ordered, solve_order{iO}];
    iO = iO + 1;
end
end

function order = recursiveSolve(opt_dep, iO)
% Recursive function to determine the solve order of the optimisation 
% problems
if ~any(opt_dep(iO,:))
    % No dependency, start of order
    order = iO;
else
    % Existing dependency
    order = [];
    ind_opt = find(opt_dep(iO,:)==1);
    for i = 1:length(ind_opt) % If length(ind_opt)>1 there are brances
        tmp_order = recursiveSolve(opt_dep, ind_opt(i));
        order = [tmp_order, order];
    end
    % Add itself to the order
    order = [order, iO];
end

end

