function [junk, M] = EigenGaps(eigenvals)
    num = size(eigenvals,1);
    eigengaps = zeros(num-1,1);
    for i = 1:num-1
        eigengaps(i) = eigenvals(i+1)-eigenvals(i);
    end
    
    figure
    plot(1:num-1,eigengaps,'.','MarkerSize',30)
    title('EigenGaps')
    grid
    [junk, M] = max(eigengaps);
end