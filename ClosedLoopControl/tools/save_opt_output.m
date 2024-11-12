function save_opt_output(x,fval,path)
%SAVE_OPT_OUTPUT saves the x value and cost function value to a csv in the
%simulation folder
writematrix(x,[path,'x.csv'])
end

