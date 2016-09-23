%p_k newton 
function p_k_newton =  P_k_newton(x1,x2)
    p_k_newton  = (-1*(inv(Hess_Rosenbrock(x1,x2)))*Grad_Rosenbrock(x1,x2));
end