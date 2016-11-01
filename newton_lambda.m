% Input: delta
% Output: lambda - the calculated value based on NW Alg 4.3
function lambda = newton_lambda(delta, B, g, lambda_0)
%for the first iteration, start with lambda = 0
lambda = lambda_0;
% NW - Alg 4.3 - Safeguard of a small number of iterations as the solution
% should  converge quickly
    for k = 1:5
        R = chol(B+lambda*eye(2)); %conduct cholesky factorization
        p_k = (B+lambda*eye(2))\(-g);
        q_k = R'\p_k;
        lambda = lambda+(norm(p_k)/norm(q_k))^2*(norm(p_k)-delta)/(delta);
        % display (lambda) % test code
    end
end

 