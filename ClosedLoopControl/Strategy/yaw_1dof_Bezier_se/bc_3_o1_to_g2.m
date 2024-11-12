function g2 = bc_3_o1_to_g2(o1, g1, rn)
%O1_TO_O2 Translation of a opimization value [0,1] into g2 (needed for bc
%and bc_d)
arguments
    o1 (:,1) double
    g1 (:,1) double
    rn (1,1) double
end

g2 = (1-o1).* ( 1/3*( -sqrt( rn^2 + 3*rn * g1) - rn + 3 * g1)) + ...
        o1 .* (-1/3*( -sqrt( rn^2 - 3*rn * g1) - rn - 3 * g1));

if sum(~isreal(g2))>0
    error('BC_3 error: Generated g2 value(s) are imaginary.')
end

end