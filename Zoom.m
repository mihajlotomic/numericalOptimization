function alpha_star = Zoom(alpha_low, alpha_high, ...
                     x, c_1, c_2, p_k)
    MAX_ITERS = 25;
    k = 0;
    while true
        alpha_j = Interpolation(alpha_high, alpha_low,...
                                x, c_1, p_k);
            
        %Compute phi(alpha_j)
        phi_alpha_j = Phi(alpha_j,x ,p_k);
        phi_alpha_low = Phi(alpha_low,x ,p_k);
        
        if ~(is_armijo_met(alpha_j,x,c_1,p_k) ||  ...
                phi_alpha_j >= phi_alpha_low)
            alpha_high = alpha_j;
        else           
            phi_prime_alpha_j = Phi_prime(alpha_j, x ,p_k);
                       
            if abs(phi_prime_alpha_j) <=  ...
                    (-1)*c_2* Phi_prime(0, x ,p_k)
            
                alpha_star = alpha_j;
                return;
            end   
        
            if phi_prime_alpha_j*(alpha_high-alpha_low) >= 0
                alpha_high = alpha_low;
            end
        
            alpha_low = alpha_j;
        end
        k = k+1;
        % Failsafe
        if k > MAX_ITERS
            alpha_star = alpha_j;
            break;
        end
    end
     
 end
