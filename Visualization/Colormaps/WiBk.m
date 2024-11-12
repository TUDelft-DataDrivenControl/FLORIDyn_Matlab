function map = WiBk(N)
% Divergent dark Red-White- dark Blue colormap 
% Data source:
% http://www.ncl.ucar.edu/Document/Graphics/ColorTables/MPL_RdBu.shtml
% Code adapted from Stephen Cobeldick
% Info on colormaps:
%   https://matplotlib.org/stable/tutorials/colors/colormaps.html
%% Example
% [X,Y,Z] = peaks(30);
% surfc(X,Y,Z)
% axis([-3,3,-3,3,-10,5])
% colormap(RdBu)
% caxis([-4,4])
% c = colorbar;
% c.Limits = [-5,2];
% c.Label.String = 'RbBu';
% 
%% Input and Output Arguments %%
%%% Inputs (*=default):
% N = NumericScalar, N>=0, an integer to define the colormap length.
%   = *[], use the length of the current figure's colormap (see COLORMAP).
%
%%% Outputs:
% map = NumericMatrix, size Nx3, a colormap of RGB values between 0 and 1.

if nargin<1 || isnumeric(N)&&isequal(N,[])
	N = size(get(gcf,'colormap'),1);
else
	assert(isscalar(N)&&isfinite(N)&&isreal(N),...
		'SC:viridis:N:NotRealFiniteScalarNumeric',...
		'First argument <N> must be a real finite numeric scalar.')
end
%% 
% r g b
raw = [1.0000    1.0000    1.0000
    0.9412    0.9412    0.9412
    0.8510    0.8510    0.8510
    0.7412    0.7412    0.7412
    0.5882    0.5882    0.5882
    0.4510    0.4510    0.4510
    0.3216    0.3216    0.3216
    0.1451    0.1451    0.1451
         0         0         0];

%
num = size(raw,1);
% With small extrapolation when N>num:
% vec = linspace(0,num+1,N+2);
% map = interp1(1:num,raw,vec(2:N+1),'linear','extrap');
vec = linspace(1,num,N);
map = interp1(1:num,raw,vec(1:N),'linear','extrap');
% Interpolation only for all values of N:
%map = interp1(1:num,raw,linspace(1,num,N),'spline')
% Range limits:
map = max(0,min(1,map));

end