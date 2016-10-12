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
        armijo_check = (Phi(alpha_k, x(1),x(2),p_k(1),p_k(2)) <= ...               
                        Phi(0,x(1),x(2),0,0) + ...
                        c_1*alpha_k*...
                        Grad_Rosenbrock(x(1),x(2))'*p_k(x(1),x(2)));
   else
       % Use the steepestdescent version
   end
end