function x = DirectLanczos(A, b, x0, maxiter, tol)
    x = x0;
    r = b-A*x;
    beta = norm(r,2); % fixed beta
    v_old = r./beta; % v1
    z = beta; % z_1(1)
    alpha = (A*v_old)'*v_old; % alpha1
    eta_old = alpha; % eta1
    eta_new = eta_old; % eta1
    w = A*v_old - alpha.*v_old; % w = Av1-alpha*v1
    beta_new = norm(w,2); % beta2 = ||w||
    v_new = w./beta_new; % v2 = w/beta2
    p = v_old./eta_new; % p = v1/eta1

    m = 1;
    disp(beta_new * abs(z/eta_new))
    while beta_new * abs(z/eta_new) > tol && m < maxiter
        m = m+1;
        w = A*v_new - beta_new.*v_old; % w = Av2 - beta2*v1
        alpha = w'*v_new; % alpha = w^T v2
        w = w - alpha.*v_new; % w = w - alpha*v2
        beta_old = beta_new; % beta_old = beta2
        beta_new = norm(w,2); % beta3 = ||w||
        v_old = v_new; % v_old = v2
        v_new = w./beta_new; % v3 = w / beta3
        lambda = beta_old/eta_old; % lambda2 = beta2 / eta1
        eta_old = eta_new; % eta_old = eta1
        eta_new = alpha - lambda*beta_old; % eta2 = alpha - lambda*beta2
        z = (-1)*lambda*z;
        p = (v_old - beta_old.*p)./eta_new;
        x = x + z.*p;
    end
end