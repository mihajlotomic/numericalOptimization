function phi_prime = Phi_prime(alpha_k, x, p_k)
   % Comments?
   phi_prime = Grad_Rosenbrock(x(1) + alpha_k*p_k(1), x(2)+alpha_k*p_k(2)'*p_k);
   
end