function minCut = randomContraction(graph)
    minCut = size(graph,1);
    N = numel(unique(graph));
    % repeatedly do random contraction
    for rep = 1:(N^2*log(N))
        trialGraph = graph;
        % keep contracting until two vertices are left
        while numel(unique(trialGraph)) > 2
            % determine two vertices
            nRow = size(trialGraph,1);
            merge = randi([1 nRow]);
            [x, y] = find(trialGraph == trialGraph(merge,2));
            
            % contract them
            trialGraph(x(y==1),1) = trialGraph(merge,1);
            trialGraph(x(y==2),2) = trialGraph(merge,1);
            
            % remove self-loops
            trialGraph(trialGraph(:,1) == trialGraph(:,2),:) = [];
        end
        
        % find the minimal cut so far
        trialMin = size(trialGraph,1);
        minCut = min([minCut, trialMin]);
        % disp([num2str(rep),': ',num2str(minCut),' (',num2str(trialMin),')'])
    end
end