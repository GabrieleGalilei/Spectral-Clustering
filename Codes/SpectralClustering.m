function [L, conncomp, eigenvals, U, M, IDX] = SpectralClustering(k, sigma, num, ...
    tol, maxit, dataset, normlapl, clustmeth, solvels)

    switch dataset
        case 'circle'
            X = load('Circle.mat');
            X = X.X;
        case 'spiral'
            X = load('Spiral.mat');
            X = X.X;
        case 'landmines'
            X = load('mixoutALL_shifted.mat');
            X = X.mixout;
            A = X{1,1};
            for i = 1:fix(length(X)/200)
                A = horzcat(A,X{1,i});
                A(:,all(~A,1)) = [];
            end
            X = A';
    end
    n = size(X,1); % n = n. of points
    
    %% 1. ADJACENCY MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dist = DistMatrix(X); 
    % it creates the distance matrix between the points
    W = WeigthMatrix(dist, sigma, k); 
    % it creates the adjecency matrix for
    % the points in the dataset, where if one point is neighbour for
    % another one the weigth is the simlarity, else is 0.
    
    G = graph(W);
    figure
    plot(G,'layout','force')
    title('Adjacency graph')
    subtitle(['k = ',num2str(k)])
    
    %% 2. DEGREE MATRIX, LAPLACIAN MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D = sparse(diag(W*ones(n,1)));
    % compute the Degree matrix
    
    switch normlapl
        case 'unnorm'
            L = sparse(D-W); % compute the Laplacian matrix
        case 'symnorm'
            Dhalf = D^(-0.5);
            L = speye(n) - Dhalf*W*Dhalf;
            L = sparse(L);
    end
    
    %% 3+4. EIGENVALUES COMPUTATION, N.o OF CLUSTERS, CONNECTED COMPONENTS
    [eigenvals, eigenvecs, conncomp] = kSmallestEigs(L, num, tol, maxit, solvels);
    % this function computes the 'num' smallest eigenvalues and the
    % corresponding eigenvectors of the Laplacian matrix; in addition, it 
    % finds the number of connected components of the adjacency graph.
    
    [~,M] = EigenGaps(eigenvals);
    % I decided to compute a suitable number of clusters using the 
    % eigengaps method.
    
    %% 5. CLUSTERS' EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U = eigenvecs(:,1:M);
    % here I recollect the eigenvectors corresponding to the M smallest 
    % eigenvalues.
    
    %% 6+7. CLUSTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch clustmeth
        case 'kmeans'
            IDX = kmeans(U,M,'Replicates',1000);
        case 'kmedoids'
            IDX = kmedoids(U,M,'Replicates',1000);
        case 'ward'
            IDX = clusterdata(U,'Linkage','ward','SaveMemory','on','Maxclust',M);
        case 'euclidean'
             IDX = clusterdata(U,'Maxclust',M);
    end
    %% 8. CLUSTER PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PlotClusters(X, M, IDX, 30, k)
    
    %[idx, V, D] = spectralcluster(X, M, 'LaplacianNormalization', 'none', ...
    %    'SimilarityGraph', 'knn', 'NumNeighbors', k, 'KNNGraphType', ...
    %    'complete', 'ClusterMethod', 'kmeans');
    %c = eigs(L);
    %m = min([num size(D,1)]);
    %disp([c(1:m) eigenvals(1:m) D(1:m)])
    
    %PlotClusters(X, m, idx, 30, k)
end