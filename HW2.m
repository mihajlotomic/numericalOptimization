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
% Steepest descent method
use_newton = 0;
p_k = P_k(x(1),x(2));

Interpolation(alpha_0,x,c_1,p_k,use_newton)





