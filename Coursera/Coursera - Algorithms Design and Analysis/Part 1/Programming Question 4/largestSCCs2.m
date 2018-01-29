function sizes = largestSCCs2(graph)
% output the sizes of the strongly connected components in a
% directed graph, each row specifying tail and head of an edge

% build search stack

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
                
                % implement DFS search as a stack
                ex(i) = 1;
                leader(i) = s;
                searchStack = i;
                searchN = 1; % search array size
                while searchN > 0
                    currentNode = searchStack(searchN); % current node in DFS
                    nextNodes = G(G(:,1) == currentNode,2); % all nodes having an edge from current
                    nextNodes = nextNodes(~ex(nextNodes)); % among which find the nodes not explored
                    
                    searchContent = numel(searchStack);
                    if ~isempty(nextNodes) % if any node is found
                        % mark it as explored and set its leader as source
                        ex(nextNodes(1)) = 1;
                        leader(nextNodes(1)) = s;
                        
                        % push it onto searchArray
                        searchN = searchN + 1;
                        if searchN > searchContent
                            copyStack = zeros(1,2*searchContent);
                            copyStack(1:searchN-1) = searchStack;
                            searchStack = copyStack;
                        end
                        searchStack(searchN) = nextNodes(1);
                        
                    else % if not found
                        % store finishing time to current node
                        t = t + 1;
                        f(currentNode)  = t;
                        searchN = searchN - 1;
                        if searchN < searchContent / 4 && searchContent >= 4
                            searchStack = searchStack(1:(searchContent/2));
                        end
                    end
                    
                end
                
            end
        end
    end

end