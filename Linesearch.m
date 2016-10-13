function alpha_star = Linesearch (alpha_0, alpha_1, x, c_1, c_2, p_k, use_newton)
    % Inputs:
    %       alpha_0, alpha_1 - 
    %       x
    %       c1, c2 - c1: close to 0, c2: close to 1
    %       use_newton - apply newton method (1) or steepest descent (0)

    % Initial conditions
    alpha_i_minus_1 = alpha_0;
    alpha_i = alpha_1;
    alpha_max = alpha_1;
    i=1;
 
    while true
        phi_alpha_i = Phi(alpha_i, x, p_k);
        if (phi_alpha_i > ( Rosenbrock(x(1),x(2)) + ....
                           c_1*alpha_i*Phi_prime(0,x,p_k)) ) || ...
           (phi_alpha_i >= Phi(alpha_i_minus_1, x, p_k) && i > 1)
            alpha_star = Zoom(alpha_i_minus_1, alpha_i,x, c_1, c_2, ...
                              p_k, use_newton);
            return;
        end
        
        phi_prime_alpha_i = Phi_prime(alpha_i, x, p_k);
        if abs(phi_prime_alpha_i) <= (-1)*c_2*Phi_prime(0,x,p_k)
            alpha_star = alpha_i;
            return;
            
        end    
        
        if phi_prime_alpha_i >=0
            alpha_star = Zoom(alpha_i, alpha_i_minus_1, x, c_1, c_2, ...
                              p_k, use_newton);
            return;
        end
        
        % Update alpha_i_minus_1 and alpha_i 
        % and resume linesearch 
        alpha_i_minus_1 = alpha_i;
        alpha_i = (alpha_max-alpha_i)/2;
        i=i+1;
   end    
           
    
end    