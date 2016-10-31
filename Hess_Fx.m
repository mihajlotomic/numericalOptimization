% gradient of HW3 function
function hess_f_x = Hess_Fx(x1, x2)
    % Output:scalar
  
    hess_f_x = [2-40*x2+120*x1.^2,-40*x1;
                    -40*x1,repmat(20, size(x1))];
             
end