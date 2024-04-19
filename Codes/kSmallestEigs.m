function [eigenvals, eigenvecs, conncomp] = kSmallestEigs(L, num, tol, maxit, solvels)
    n = size(L,1);
    eigenvals = zeros(num,1);
    eigenvecs = zeros(n,num);
    A = L;
    found = 0;
    
    for k = 1:num
        [lambda, v] = inversepower(A, tol, maxit, solvels);
        eigenvals(k) =  lambda;
        A = deflation(A, v);
        eigenvecs(k:n,k) = v;
        if lambda > 1.0e-16 && found == 0
            found = k;
        end
    end
    
    [eigenvals,I] = sort(eigenvals);
    eigenvecs = eigenvecs(:,I);
    
    conncomp = found-1;
end