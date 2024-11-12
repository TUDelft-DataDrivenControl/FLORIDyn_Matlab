function bc_d = bc_3_d(t_n,g0,g1,g2,g3)
%BC_3_D is the derivative of a cubic bezier-curve, with t_n \in [0,1].
%   Returns a matrix of the size (timesteps x turbines)
arguments
    t_n (:,1) double
    g0 (1,:) double
    g1 (1,:) double
    g2 (1,:) double
    g3 (1,:) double
end

bc_d = 3*(1-t_n).^2  * (g1 - g0)  +...
    6*(1-t_n).* t_n * (g2 - g1) +...
    3*t_n.^2 * (g3 - g2);

end

