function p_k_samp = Samp_P_k(x1,x2)
%    y = (x1+x2^2)^2
%    grad_out = [2*(x1+x2^2);
%                4*x2*(x1+x2^2)];
% Steepest Descent 
p_k_samp = (-1*(grad_f_x(x1,x2))/norm(grad_f_x(x1,x2)));

end