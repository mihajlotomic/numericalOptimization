function phi= Phi(alpha_k, x, p_k)
    %Inputs: alpha_k - step length (scalar)
    %        x: 2x1 
    %        p_k: 2x1
    % Output: Scalar - f(x+alpha*p_k)
    phi= Rosenbrock(x(1) +alpha_k*p_k(1), x(2)+alpha_k*p_k(2));
   
end