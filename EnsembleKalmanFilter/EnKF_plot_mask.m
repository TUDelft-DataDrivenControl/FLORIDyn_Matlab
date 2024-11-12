function EnKF_plot_mask(X,Y,M,levels, color)
%ENKF_PLOT_MASK Plots a mask across a given flow field based on the values
%of M
%   Elements of M need to be between 0 and 1, 1 meaning they should be
%   trused the most.
% ======================================================================= %
% Example Code:
%
% [X,Y] = meshgrid(1:20,1:20);
% M = rand(20,20);
% figure
% contourf(X,Y,M,20,'EdgeColor','none')
% hold on
% EnKF_plot_mask(X,Y,M,.8, [1,1,1])
% hold off
% ======================================================================= %

arguments
    X (:,:) double % X Coordinates of map
    Y (:,:) double % Y Coordinates of map
    M (:,:) double % "Truth" values of map
    levels (1,:) double % Levels
    color (1,3) double = 0 % Color of map
end

%% Normalise M & levels to range between 0 and 1
% Normalise levels to the range of M from low to high.
levels  = sort((levels-min(M,[],"all"))/(max(M,[],"all")-min(M,[],"all")));

% Normalise M
M       = (M-min(M,[],"all"))/(max(M,[],"all")-min(M,[],"all"));

%% Plot

% Define the white color with varying opacity
if length(levels)>2
    face_alpha = levels; % Opacity for each level
elseif length(levels)==2
    face_alpha = [0,1];  % Opacity for each level
else
    face_alpha = 1;
end

% Set color and opacity for each contour level
for i = 1:length(levels)
    contourf(X, Y, M, [levels(i), levels(i)], 'FaceAlpha',face_alpha(i),...
        'FaceColor', color,'EdgeColor','none');
end

end

