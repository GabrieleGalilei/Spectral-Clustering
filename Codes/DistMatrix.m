function dist = DistMatrix(X)
    [n,dim] = size(X);
    dist = zeros(n,n);
    
    for i = 1:n 
        for j = i+1:n 
            for d = 1:dim
                dist(i,j) = dist(i,j) + (X(i,d)-X(j,d))^2;
            end
        end
    end
    
    dist = dist+dist';
end