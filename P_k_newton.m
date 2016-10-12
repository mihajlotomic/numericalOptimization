%p_k newton 
function p_k_newton =  P_k_newton(x1,x2)
    % 2x2 times 2x1 returns 2x1
    % b = Ax -> inv(A)*b = A\b = x
    % ToDo: Use A\b in future version 
    %p_k_newton  = (-1*(inv(Hess_Rosenbrock(x1,x2)))*Grad_Rosenbrock(x1,x2));
    p_k_newton  = -1*(Hess_Rosenbrock(x1,x2)) \ Grad_Rosenbrock(x1,x2);
end