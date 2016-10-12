%% Steepest Descent
clear
MAXITER = 10^4;
alpha = 1;
x=[-1.2;1];
c=10^(-4);
rho=1/2;
Grad_Rosenbrock(x(1),x(2));

k = 0;
display(sprintf('x^T \t\t\t\t\tRosenbrock \t\t\t p^SD \t\t\t\t\talpha \t\titer'));

%stopping condition
while (abs(Rosenbrock(x(1),x(2)))) > (10^(-8))
%search direction
    p_k = P_k(x(1),x(2));
    
    %check armijo condition
    while (Rosenbrock(x(1)+alpha*p_k(1),x(2)+alpha*p_k(2))) > ...
         ( Rosenbrock(x(1),x(2)) + ... 
           c*alpha*P_k(x(1),x(2))'*Grad_Rosenbrock(x(1),x(2))  ...
         )   
     
     % update alpha
        alpha=rho*alpha;
    end
    
    %update the minimizer
    x = x + alpha*P_k(x(1),x(2));
   alpha = 1;
   k=k+1;
   if mod(k,1000) == 0
           display(sprintf('[%1.6f,%1.6f]\t\t[%1.6f]\t\t[%1.6f,%1.6f]\t\t%1.6f\t%1.6f',...
       x(1),x(2),Rosenbrock(x(1),x(2)), P_k(x(1),x(2)) , alpha, k));
   end
   if k > MAXITER
       break;
   end
end
display(sprintf('[%1.6f,%1.6f]\t\t[%1.6f]\t\t[%1.6f,%1.6f]\t\t\t%1.6f\t%1.6f',...
       x(1),x(2),Rosenbrock(x(1),x(2)), P_k(x(1),x(2)) , alpha, k));
%% Newton Method Solution
%
clear
alpha = 1;
x=[-1.2;1];
c=10^(-4);
rho=1/2;
Grad_Rosenbrock(x(1),x(2));

k = 0;
display(sprintf('x^T \t\t\t\t\t\t\tRosenbrock \t\t\t\t p^N \t\t\t\t\t\t\t\t alpha \t\t\titer'));

while (norm(Grad_Rosenbrock(x(1),x(2)))) > (10^(-8))
    k=k+1;

    p_k_newton = P_k_newton(x(1),x(2));
    
    
    while (Rosenbrock(x(1)+alpha*p_k_newton(1),x(2)+alpha*p_k_newton(2))) > ...
         ( Rosenbrock(x(1),x(2)) + ... 
           c*alpha*P_k_newton(x(1),x(2))'*Grad_Rosenbrock(x(1),x(2))  ...
         )  

        alpha=rho*alpha;
    end
   
    x = x + alpha*P_k_newton(x(1),x(2));
    display(sprintf('[%1.6f,%1.6f]\t\t\t[%1.6f]\t\t\t\t[%1.6f,%1.6f]\t\t\t\t\t%1.6f\t\t%1.6f',...
         x(1),x(2),Rosenbrock(x(1),x(2)), P_k_newton(x(1),x(2)) , alpha, k));
    alpha=1;
end

       

   


