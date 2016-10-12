% Steepest descent direction via gradient of rosenbrock
function grad_rosenbrock = Grad_Rosenbrock(x1, x2)
    % Output: 2x1 partial derivatives
  
    grad_rosenbrock = [-400*(x2-x1^2)*(x1)-2*(1-x1) ;
                        200*(x2-x1^2) ];
end