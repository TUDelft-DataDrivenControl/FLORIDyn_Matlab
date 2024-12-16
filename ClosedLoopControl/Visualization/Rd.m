function map = Rd(N)
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
raw = [0.4155    0.0037    0.1234
    0.4385    0.0111    0.1271
    0.4616    0.0185    0.1308
    0.4847    0.0258    0.1345
    0.5077    0.0332    0.1382
    0.5308    0.0406    0.1419
    0.5539    0.0480    0.1456
    0.5769    0.0554    0.1492
    0.6000    0.0627    0.1529
    0.6115    0.0664    0.1548
    0.6461    0.0775    0.1603
    0.6577    0.0812    0.1622
    0.6923    0.0923    0.1677
    0.7008    0.0997    0.1712
    0.7174    0.1329    0.1869
    0.7230    0.1439    0.1922
    0.7396    0.1772    0.2078
    0.7506    0.1993    0.2183
    0.7562    0.2104    0.2235
    0.7728    0.2436    0.2392
    0.7839    0.2657    0.2497
    0.7949    0.2879    0.2601
    0.8005    0.2990    0.2654
    0.8171    0.3322    0.2810
    0.8281    0.3543    0.2915
    0.8392    0.3765    0.3020
    0.8438    0.3871    0.3101
    0.8577    0.4189    0.3346
    0.8669    0.4401    0.3509
    0.8761    0.4614    0.3672
    0.8807    0.4720    0.3753
    0.8946    0.5038    0.3998
    0.9038    0.5250    0.4161
    0.9130    0.5463    0.4324
    0.9223    0.5675    0.4487
    0.9315    0.5887    0.4650
    0.9407    0.6099    0.4813
    0.9453    0.6205    0.4894
    0.9576    0.6512    0.5151
    0.9603    0.6678    0.5363
    0.9631    0.6844    0.5576
    0.9659    0.7010    0.5788
    0.9686    0.7176    0.6000
    0.9714    0.7343    0.6212
    0.9742    0.7509    0.6424
    0.9755    0.7592    0.6531
    0.9797    0.7841    0.6849
    0.9825    0.8007    0.7061
    0.9852    0.8173    0.7273
    0.9880    0.8339    0.7486
    0.9908    0.8505    0.7698
    0.9912    0.8631    0.7878
    0.9894    0.8717    0.8025
    0.9885    0.8760    0.8099
    0.9857    0.8890    0.8321
    0.9839    0.8976    0.8468
    0.9820    0.9062    0.8616
    0.9802    0.9148    0.8764
    0.9783    0.9234    0.8911
    0.9765    0.9320    0.9059
    0.9746    0.9406    0.9206
    0.9737    0.9449    0.9280
    0.9709    0.9579    0.9502
    0.9691    0.9665    0.9649];

%
num = size(raw,1);
% With small extrapolation when N>num:
vec = linspace(1,num,N);
map = interp1(1:num,raw,vec(1:N),'linear','extrap');
% vec = linspace(0,num+1,N+2);
% map = interp1(1:num,raw,vec(2:N+1),'linear','extrap');
% Interpolation only for all values of N:
%map = interp1(1:num,raw,linspace(1,num,N),'spline')
% Range limits:
map = max(0,min(1,map));

end