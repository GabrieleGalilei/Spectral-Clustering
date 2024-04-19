function x = gmres_sparse(A, b, restart, max_iter, tol)
    n = size(A, 1);
    m = min(restart, n);

    x = zeros(n, 1);
    res = norm(b - A*x);
    V = zeros(n, m + 1);
    H = zeros(m + 1, m);

    iter = 0;
    while res > tol && iter < max_iter
        iter = iter + 1;

        V(:, 1) = (b - A * x) / res;
        s = zeros(m + 1, 1);
        s(1) = res;

        for j = 1:m
            w = A * V(:, j);
            for i = 1:j
                H(i, j) = w' * V(:, i);
                w = w - H(i, j) * V(:, i);
            end
            H(j + 1, j) = norm(w);
            if H(j + 1, j) < 1e-12
                break;
            end
            V(:, j + 1) = w / H(j + 1, j);
            % Apply Givens rotation
            for i = 1:j - 1
                temp = c(i) * H(i, j) + s(i) * H(i + 1, j);
                H(i + 1, j) = -s(i) * H(i, j) + c(i) * H(i + 1, j);
                H(i, j) = temp;
            end
            [c(j), s(j)] = givens_rotation(H(j, j), H(j + 1, j));
            temp = c(j) * H(j, j) + s(j) * H(j + 1, j);
            H(j + 1, j) = -s(j) * H(j, j) + c(j) * H(j + 1, j);
            H(j, j) = temp;
            s(j + 1) = -s(j) * H(j + 1, j);
            H(j + 1, j) = c(j) * H(j + 1, j);
        end

        y = H(1:j, 1:j) \ s(1:j);
        x = x + V(:, 1:j) * y;
        res = norm(b - A * x);
    end
end