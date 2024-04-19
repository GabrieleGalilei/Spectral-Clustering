function [A,P] = deflation(A, v)
    [m,n] = size(A);
    if m ~= n
        disp('The matrix is not square');
        return;
    end
    disp('computing deflation')
    v = v/norm(v,2);
    u = v;
    u(1) = u(1) + sign(v(1));

    w = u';
    w = (u*w)/(norm(u,2)^2);
    P = speye(n) - 2*w;
    B = P*A*transpose(P);
    A = B(2:n,2:n);
end