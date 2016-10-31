% gradient of HW3 function
function grad_f_x = Grad_Fx(x1, x2)
    % Output: 2x1 partial derivatives
  
    grad_f_x = [-40.*(x2-x1.^2).*(x1)-2.*(1-x1);
                        20.*(x2-x1.^2) ];
end