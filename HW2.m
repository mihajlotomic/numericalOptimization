%% 
clear
TOLERANCE = 10^-8;
MAXITER = 20e3;
alpha_val = 1; % alpha_0, alpha_max
x=[-1.2;1.2];
c1=10^(-4);
c2= 9/10;

Grad_Rosenbrock(x(1),x(2));

k = 0;
fprintf('x^T \t\t\t\t\tRosenbrock \t\t\t p^ \t\t\t\t\talpha \t\titer\n');

while abs(Rosenbrock(x(1),x(2)))  > TOLERANCE 
    % Search direction
    %p_k = P_k(x(1),x(2));
    p_k = P_k_newton(x(1), x(2));

    alpha_val = Linesearch(alpha_val, x, c1, c2, p_k);
    
    % Update the minimizer
    x = x + alpha_val*P_k_newton(x(1),x(2));

    k=k+1;
    if mod(k,1) == 0
        fprintf('[%1.6f,%1.6f]\t\t[%1.6f]\t\t[%1.6f,%1.6f]\t\t%1.6f\t%1.6f\n',...
            x(1),x(2),Rosenbrock(x(1),x(2)), P_k_newton(x(1),x(2)) , alpha_val, k);
    end
    alpha_val = 1; % For Newton method
    if k > MAXITER
        error('Maximum iteration met - abort search', MAXITER)
    end
end
fprintf('[%1.6f,%1.6f]\t\t[%1.6f]\t\t[%1.6f,%1.6f]\t\t\t%1.6f\t%1.6f\n',...
       x(1),x(2),Rosenbrock(x(1),x(2)), P_k_newton(x(1),x(2)) , alpha_val, k);