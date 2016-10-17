function alpha_star  = Interpolation(alpha_k, alpha_low, ...
                                       x, c_1, p_k)
    epsilon_1 = 1e-3;
    epsilon_2 = 1e-3;
    % Persistent as a way to keep track of the number of function calls
    persistent j; 

    % Step length is bound by alpha_low
    if alpha_k <= alpha_low
        alpha_star = alpha_low;
        return;
    end
    
    if is_armijo_met(alpha_k,x,c_1,p_k)
        alpha_star = alpha_k;
        % Terminate serach - we've found our step length
        return;
    end
    % calculate alpha1
    alpha_1 =  -1*Phi_prime(0, x, p_k)*(alpha_k^2) / ...
               (2*(Phi(alpha_k, x ,p_k) - Phi(0, x, p_k) - ...
                  Phi_prime(0, x, p_k)*alpha_k));

    if is_armijo_met(alpha_1,x,c_1,p_k)
        alpha_star = alpha_1;
        % Terminate serach - we've found our step length
        return;
    end
    
   %breaking [a;b] matrix into 3 pieces: k,l,m 
    k = 1/((alpha_k^2)*(alpha_1^2)*(alpha_1-alpha_k)); 
    l =[alpha_k^2,  -1*(alpha_1^2);
        -1*(alpha_k^3),  alpha_1^3];
    m = [Phi(alpha_1, x, p_k) - Phi(0, x, p_k) - Phi_prime(0, x, p_k)*alpha_1;...
         Phi(alpha_k, x, p_k) - Phi(0, x, p_k) - Phi_prime(0, x, p_k)*alpha_k ...
        ];
        
    ab=k*l*m;
    a=ab(1);
    b=ab(2);
    % calculate alpha 2
    alpha_2 = -1*b+sqrt(b^2-3*a*Phi_prime(0, x, p_k))/(3*a);
    while ~is_armijo_met(alpha_2,x,c_1,p_k)
        % NEED A SAFEGUARD!!!   
        if abs(alpha_2 - alpha_1) < epsilon_1...
            || abs(alpha_2) < epsilon_2
              alpha_2 = alpha_1 / 2;
              
        end
        if isempty(j)
            j = 0;
        else
            j = j+1;
        end
        if j > 25
            alpha_star = alpha_1 / 2;
            j = 0;
            return;
        end
        Interpolation(alpha_2,alpha_low,x,c_1,p_k);

    end
    % Terminate serach - we've found our step length
    alpha_star = alpha_2;     
end
