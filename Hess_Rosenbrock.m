% Steepest descent direction via gradient of rosenbrock
function hess_rosenbrock = Hess_Rosenbrock(x1,x2)
hess_rosenbrock = [2-400*x2+1200*x1^2,-400*x1;
                    -400*x1,200 ];
%hess_rosenbrock = [ diff(diff(Rosenbrock,x1),x1) diff(diff(Rosenbrock,x1),x2);...
        %diff(diff(Rosenbrock,x2),x1) diff(diff(Rosenbrock,x2),x2) ];
end