%p_k steepest descent 
function p_k =  P_k(x1,x2)
    % Norm() returns a scalar, Grad_Rosenbrock a 2x1 vector, so 
    % output is 2x1
    p_k  = (-1*(Grad_Rosenbrock(x1,x2)))/norm(Grad_Rosenbrock(x1,x2));
end