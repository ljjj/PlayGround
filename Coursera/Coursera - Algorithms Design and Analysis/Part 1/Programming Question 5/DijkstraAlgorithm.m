function A = DijkstraAlgorithm(graph, source)
    % initialization
    n = max(max(graph(:,1:2))); % total number of nodes
    X = source; % vertices processed so far
    A = zeros(1,n); % computed shortest path distances
    maxPath = 1000000;
    A = A + maxPath;
    A(source) = 0;
    
    while numel(X) < n
        potentialPaths = [];
        for x = 1:numel(X)
            additionalPaths = [];
            for j = 1:2
                additionalPaths = [additionalPaths; graph(graph(:,j) == X(x),[j 3-j 3])];
            end
            additionalPaths(:,3) = additionalPaths(:,3) + A(X(x));
            potentialPaths = [potentialPaths; additionalPaths];
        end
        
        [addP, addInd] = min(potentialPaths(:,3));
        addX = potentialPaths(addInd,2);
        X = [X addX];
        A(addX) = addP;
        
        for x = 1:numel(X)
            for j = 1:2
                graph(graph(:,j) == X(x) & graph(:,3-j) == addX,:) = [];
            end
        end
    end
end