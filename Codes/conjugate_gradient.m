function x = conjugate_gradient(A, b)
% A: sparse symmetric matrix
% b: right-hand side vector

% set convergence tolerance and maximum number of iterations
tol = 1e-6;
max_iter = 1000;

% initialize solution vector and residual
n = size(A, 1);
x = zeros(n, 1);
r = b - A * x;

% initialize search direction and residual norm
p = r;
r_norm = norm(r);

% loop until convergence or maximum number of iterations reached
for iter = 1:max_iter
    % compute step size
    Ap = A * p;
    alpha = (r_norm^2) / (p' * Ap);
    
    % update solution vector and residual
    x = x + alpha * p;
    r = r - alpha * Ap;
    
    % check convergence
    r_norm_new = norm(r);
    if r_norm_new < tol
        break;
    end
    
    % update search direction
    beta = (r_norm_new^2) / (r_norm^2);
    p = r + beta * p;
    r_norm = r_norm_new;
end

end