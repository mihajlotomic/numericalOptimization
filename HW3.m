x_k = [0;-1];
k = 1;
p_k = p_k_dl(x_k(1),x_k(2));
TOLERANCE = 10^-8;
eta = 0.125;
deltak = (0:.25:2);

% 
Grad_Fx(x_k(1), x_k(2));
Fx(x_k(1),x_k(2));
Hess_Fx(x_k(1),x_k(2));

while abs(Fx(x_k(1),x_k(2)))  > TOLERANCE 
    while Fx(x_k(1)+p_k(1),x_k(2)+p_k(2)) > Fx(x_k(1),x_k(2)) % IS THIS CORRECT????
        p_k = p_k_dl(x_k(1),x_k(2));
        m_k_0 = m_k(0,0,p_k);
        m_k_p_k = m_k(x_k(1),x_k(2),p_k);

        rho_k = (Fx(x_k(1),x_k(2))-Fx(x_k(1)+pk(1),x_k(2)+pk(2)))/(m_k_0-m_k_p_k);

        if rho_k < 1/4
            deltak_plus1=1/4*deltak;
        else
            if rho_k > 3/4 && norm(p_k_dl) == deltak
            deltak_plus1 = min(2*deltak,2);
            else
                deltak_plus1 = deltak;
            end
        end  
            
        if rho_k > eta
            x_k = x_k+pk;
%%need to understand this better
        else
            x_k = x_k;
        end
        k=k+1;
    
    end
end