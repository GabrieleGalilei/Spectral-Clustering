function PlotClusters(X, M, idx, sz, k)
    [n,dim] = size(X);
    colors = cool(M);
    
    if dim == 2
        figure
        for i = 1:n
            title('Clusters')
            subtitle(['k = ',num2str(k),', n.clusters = ',num2str(M)])
            scatter(X(i,1), X(i,2), sz, colors(idx(i)), 'filled')
            hold on
        end
        grid
    elseif dim == 3
        figure
        for i = 1:n
            title('Clusters')
            subtitle(['k = ',num2str(k),', n.clusters = ',num2str(M)])
            scatter3(X(i,1), X(i,2), X(i,3), 30, colors(idx(i)), 'filled')
            hold on
        end
        grid
    end