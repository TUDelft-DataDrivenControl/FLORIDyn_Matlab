% === Ramp event wind direction ===
varphi_0    = 225;
varphi_del  = 10;
t_0         = 31000;
t_1         = 31600;
t           = linspace(t_0,t_1);

varphi = varphi_del/2*(sin(pi*(t-t_0)/(t_1-t_0) - pi/2)+1)+varphi_0;


wind_dir = ...
    [0, varphi_0
    t', varphi';
    90000, varphi_0 + varphi_del];