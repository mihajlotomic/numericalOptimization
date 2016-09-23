alpha = 1;
x=[1;1];
c=10^(-4);
rho=1/2;
p_k = Samp_P_k(x(1),x(2));

count = 0 

% Steepest Descent
display(sprintf('x^T \t\t\t\tf(x) \t\t\t p^SD \t\t\t alpha'));
while (f_x(x(1)+alpha*p_k(1),x(2)+alpha*p_k(2))) > ...
         (f_x(x(1),x(2)) + alpha*Samp_P_k(x(1),x(2))'*grad_f_x(x(1),x(2)))        

   display(sprintf('[%1.6f,%1.6f]\t\t[%1.6f]\t\t[%1.6f,%1.6f]\t\t%1.6f',x(1),x(2),f_x(x(1),x(2)), Samp_P_k(x(1),x(2)) , alpha));
   if count <= 1
       alpha = 1;
   else
       alpha = 1/2*1/(2^count);
   end
    x = x + alpha*Samp_P_k(x(1),x(2));
    p_k = P_k(x(1),x(2));
    %alpha=rho*alpha;
%   % TEST CODE
%     if ( mod(count,COUNT) == 0) 
%         sprintf('Iteration = %d',count)
%     end
      count = count + 1;
      %if f_x(x(1),x(2)) < 10^-8 
      %    break;
      %end
end