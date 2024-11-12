function l = getLayout(n)
% Generates a subplot layout for n plots
f = factor(n);
switch length(f)
    case 1
        l = [f,1];
    case 2
        l = [max(f),min(f)];
    otherwise
        l = [prod(f(1:2:end)), prod(f(2:2:end))];
        l = [max(l),min(l)];
end
end

