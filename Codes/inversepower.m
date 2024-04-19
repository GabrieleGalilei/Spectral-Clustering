function [lambda, v, k] = inversepower(A, tol,  numitr, solvels) 
    [m,n] = size(A);

    if m ~= n
        disp('matrix A is not square');
        return;
    end
    disp('computing inverse power')

    lambda = 1;
    invlambda_new = 1/lambda;
    v_old = ones(n,1)./norm(ones(n,1),2);

    for k = 1:numitr
        invlambda_old = invlambda_new;
        switch solvels
            case 'conjgrad'
                v_new = conjugate_gradient(A, v_old);
            case 'backslash'
                v_new = A\v_old;
            case 'gmres'
                v_new = gmres_sparse(A, v_old, 20, 100, 1.0e-07);
            case 'lanczos'
                v_new = DirectLanczos(A, v_old, v_old, numitr, tol);
        end
        invlambda_new = v_old'*v_new;
        v_old = v_new./norm(v_new,2);
        if abs(invlambda_new-invlambda_old) < tol
            break;
        end
    end
    lambda = 1/invlambda_new;
    v = v_old;
    return;
end