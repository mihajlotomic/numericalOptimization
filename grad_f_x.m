function grad_out = grad_f_x(x1,x2)

%    y = (x1+x2^2)^2
    grad_out = [2*(x1+x2^2);
                4*x2*(x1+x2^2)];
end
