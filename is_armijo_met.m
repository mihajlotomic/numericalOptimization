% Armijo Condition - 3.56 
% Inputs: step_length (alpha_k)
%         x (size: 2x1)
%         c_1 constant, usually 1e-4
%         p_k step direction (2x1)
%         use_newton - Boolean to use Newton method
% Output: True/False that evaluates the Armijo Condition
function armijo_check = is_armijo_met(alpha_k,x,c_1,p_k, use_newton)
   % f(x+alpha*p_k) = phi(alpha)
   % f(k)           = phi(0)
   % grad_f_k' p_k  = phi_prime(0)
   if use_newton
        armijo_check = (Phi(alpha_k, x, p_k) <= ...               
                        Phi(0, x, p_k) + c_1*alpha_k*Phi_prime(0, x, p_k));
   else
       % Use the steepestdescent version
   end
end