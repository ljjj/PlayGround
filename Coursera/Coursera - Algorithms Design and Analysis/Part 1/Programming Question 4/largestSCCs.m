function sizes = largestSCCs(graph)
% output the sizes of the strongly connected components in a
% directed graph, each row specifying tail and head of an edge

% naive recursion

    % initialization
    n = max(max(graph)); % total number of nodes
    t = 0; % number of nodes processed so far
    s = 0; % current source vertex
    ex = zeros(1,n); % is it explored?
    leader = zeros(1,n); % leader node of the SCC where i is
    f = zeros(1,n); % finishing time of the node

    % first pass to compute "magical ordering" of nodes
    DFSLoop(flip(graph,2));
    
    % create new graph
    newGraph = graph;
    for ni = 1:n
        newGraph(graph == ni) = f(ni);
    end
    ex = zeros(1,n); % remember to reset explored flag
    
    % second pass to discover SCCs
    DFSLoop(newGraph);
    
    % determine the size of SCCs
    sizes = histc(leader,unique(leader));
    
    % a single DFS pass
    function DFSLoop(G)        
        for i = n:-1:1
            if ~ex(i)
                s = i;
                DFS(G,i);                
            end
        end
    end

    function DFS(G, di)
        ex(di) = 1;
        leader(di) = s;
        nextNodes = G(G(:,1) == di,2); % all nodes having an edge from current
        for ds = 1:numel(nextNodes)
            if ~ex(nextNodes(ds))
                DFS(G, nextNodes(ds));
            end
        end
        t = t + 1;
        f(di) = t;
    end

end