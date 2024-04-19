function W = WeigthMatrix(dist, sigma, k)
    n = size(dist,1);
    W = zeros(n,n);
    for i=1:n 
        [~,I] = mink(dist(i,:), k+1); % find the k smallest element of the row
        % +1 because 0 on the diagonal has to be ignored
        for j = 2:k+1 % it ignores the smallest element (0)
            W(i,I(j)) = exp(-dist(i,I(j))/(2*sigma^2)); 
            % substitutes weights of neighbors with similarity 
            if W(I(j),i) == 0
                W(I(j),i) = W(i,I(j)); % make sure that the matrix is symmetric
            end
        end
        W(i,i) = 0; % it makes sure that elements on the diagonal are zero
        % but is not necessary
    end
    W = sparse(W);
end