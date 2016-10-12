function [alpha_star]  = Interpolation(alpha_k, alpha_low, ...
                                       x, c_1, p_k, use_newton)
    epsilon_1 = 1e-5;
    epsilon_2 = 10^(-5);
    
    % Step length is bound bu alpha_low
    if alpha_k <= alpha_low
        alpha_star = alpha_low;
        return;
    end
    
    if is_armijo_met(alpha_k,x,c_1,p_k,use_newton)
        alpha_star = alpha_k;
        % Terminate serach - we've found our step length
        return;
    end
    % calculate alpha 1
    alpha_1 = -1*(Grad_Rosenbrock(x(1),x(2))'*p_k(x(1),x(2))*alpha_0^2)/...
        2*(Rosenbrock(x(1)+alpha_0*p_k(1),x(2)+alpha_0*p_k(2)));

    if is_armijo_met(alpha_1,x,c_1,p_k,use_newton)
        alpha_star = alpha_1;
        % Terminate serach - we've found our step length
        return;
    end
    
   %breaking [a;b] matrix into 3 pieces: k,l,m 
    k = 1/(alpha_0^2 *alpha_1^2*(alpha_1-alpha_0)); 
    l =[alpha_0^2  -1*(alpha_1^2);
        -1*(alpha_0^3)  alpha_1^3];
    m = [(Rosenbrock(x(1)+alpha_1*p_k(1),x(2)+alpha_1*p_k(2))) - ...
         Rosenbrock(x(1),x(2))-Grad_Rosenbrock(x(1),x(2))'*p_k(x(1),x(2))*alpha_1;
         Rosenbrock(x(1)+alpha_0*p_k(1),x(2)+alpha_0*p_k(2)) - ...
         Rosenbrock(x(1),x(2))-Grad_Rosenbrock(x(1),x(2))'*p_k(x(1),x(2))*alpha_0];
        
    ab=k*l*m;
    a=ab(1);
    b=ab(2);
    % calculate alpha 2
    alpha_2 = -1*b+ sqrt(b^2-3*a*Grad_Rosenbrock(x(1),x(2))'*p_k(x(1),x(2)))/3*a;
    while ~is_armijo_met(alpha_2,x,c_1,p_k,use_newton)
        % NEED A SAFEGUARD!!!   
        if abs(alpha_2 - alpha_1) < epsilon_1...
            || abs(alpha_2) < epsilon_2
              alpha_2 = alpha_1 / 2;
        end
        Interpolation(alpha_2,x,c_1,p_k,use_newton);
    end
    % Terminate serach - we've found our step length
    alpha_star = alpha_2;
    % calculate alpha_k+1
%     alpha_k+1 = -1*(Grad_Rosenbrock(x(1),x(2))'*p_k(x(1),x(2))*alpha_0^2)/...
%         2*(Rosenbrock(x(1)+alpha_0*p_k(1),x(2)+alpha_0*p_k(2)));
%      
end

 




