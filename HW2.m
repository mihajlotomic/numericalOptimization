x=[-1.2;1];
alpha_0=1;
c_1=10^(-4);
c_2= 9/10;
%alpha_l = 
%alpha_h =
% Netwon method
use_newton = 1;
p_k = P_k_newton(x(1),x(2));

alpha_star = Interpolation(alpha_0,x,c_1,p_k,use_newton);

%%
%% Steepest Descent
clear
TOLERANCE = 10^-8;
MAXITER = 10^4;
alpha_vals = [0,1]; % alpha_0, alpha_max
x=[-1.2;1];
c1=10^(-4);
c2= 9/10;

% Netwon method
use_newton = 1;

Grad_Rosenbrock(x(1),x(2));

k = 0;
fprintf('x^T \t\t\t\t\tRosenbrock \t\t\t p^SD \t\t\t\t\talpha \t\titer\n');


while abs(Rosenbrock(x(1),x(2)))  > TOLERANCE 
    % Search direction
    display(k)
    p_k = P_k(x(1),x(2));

    Linesearch(alpha_vals(1), alpha_vals(2), x, c1, c2, p_k, 1);
    
   % Update the minimizer
   x = x + alpha_vals*P_k(x(1),x(2));
   k=k+1;
   if mod(k,1000) == 0
           fprintf('[%1.6f,%1.6f]\t\t[%1.6f]\t\t[%1.6f,%1.6f]\t\t%1.6f\t%1.6f',...
       x(1),x(2),Rosenbrock(x(1),x(2)), P_k(x(1),x(2)) , alpha, k);
   end
   if k > MAXITER
       break;
   end
end

% 
fprintf('[%1.6f,%1.6f]\t\t[%1.6f]\t\t[%1.6f,%1.6f]\t\t\t%1.6f\t%1.6f',...
       x(1),x(2),Rosenbrock(x(1),x(2)), P_k(x(1),x(2)) , alpha, k);