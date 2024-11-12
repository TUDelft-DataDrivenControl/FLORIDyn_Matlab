function map = RdBu(N)
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
raw = [0.9599    0.9645    0.9670
    0.9482    0.9589    0.9649
    0.9365    0.9534    0.9627
    0.9248    0.9479    0.9606
    0.9131    0.9423    0.9584
    0.9014    0.9368    0.9562
    0.8897    0.9313    0.9541
    0.8780    0.9257    0.9519
    0.8664    0.9202    0.9498
    0.8547    0.9146    0.9476
    0.8430    0.9091    0.9455
    0.8371    0.9063    0.9444
    0.8196    0.8980    0.9412
    0.8002    0.8882    0.9356
    0.7809    0.8784    0.9301
    0.7615    0.8685    0.9246
    0.7421    0.8587    0.9190
    0.7227    0.8488    0.9135
    0.7033    0.8390    0.9080
    0.6840    0.8291    0.9024
    0.6646    0.8193    0.8969
    0.6452    0.8095    0.8913
    0.6258    0.7996    0.8858
    0.6065    0.7898    0.8803
    0.5871    0.7799    0.8747
    0.5665    0.7687    0.8685
    0.5422    0.7533    0.8602
    0.5300    0.7456    0.8561
    0.4936    0.7226    0.8436
    0.4693    0.7072    0.8353
    0.4450    0.6918    0.8270
    0.4207    0.6764    0.8187
    0.3964    0.6611    0.8104
    0.3721    0.6457    0.8021
    0.3478    0.6303    0.7938
    0.3235    0.6149    0.7855
    0.2992    0.5995    0.7772
    0.2749    0.5842    0.7689
    0.2575    0.5696    0.7612
    0.2471    0.5557    0.7541
    0.2366    0.5419    0.7470
    0.2261    0.5280    0.7399
    0.2157    0.5142    0.7329
    0.2105    0.5073    0.7293
    0.1948    0.4865    0.7187
    0.1843    0.4727    0.7116
    0.1739    0.4588    0.7046
    0.1634    0.4450    0.6975
    0.1529    0.4311    0.6904
    0.1425    0.4173    0.6834
    0.1320    0.4035    0.6763
    0.1230    0.3875    0.6572
    0.1143    0.3709    0.6341
    0.1057    0.3543    0.6111
    0.0971    0.3377    0.5880
    0.0885    0.3211    0.5649
    0.0799    0.3045    0.5419
    0.0713    0.2879    0.5188
    0.0627    0.2713    0.4957
    0.0584    0.2630    0.4842
    0.0454    0.2381    0.4496
    0.0368    0.2215    0.4265
    0.0282    0.2048    0.4035
    0.0196    0.1882    0.3804];

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