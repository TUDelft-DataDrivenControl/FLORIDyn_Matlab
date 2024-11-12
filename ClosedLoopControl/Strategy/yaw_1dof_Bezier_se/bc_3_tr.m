function tr = bc_3_tr(t_n,g0,g1,g2,g3)
%BC_3 is the trajectory of a cubic bezier-curve, with t_n \in [0,1].
%   Returns a matrix of the size (timesteps x turbines)
arguments
    t_n (:,1) double
    g0 (1,:) double
    g1 (1,:) double
    g2 (1,:) double
    g3 (1,:) double
end

tr = (1-t_n).^3*g0 + 3*(1-t_n).^2 .* t_n * g1 + ...
    3*(1-t_n).* t_n.^2 * g2 + t_n.^3 * g3;

end

