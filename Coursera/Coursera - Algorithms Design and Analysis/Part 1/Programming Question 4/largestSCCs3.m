function sizes = largestSCCs3(graph)
% output the sizes of the strongly connected components in a
% directed graph, each row specifying tail and head of an edge

% build search array with path stack to avoid looking for the same next
% nodes repeatedly

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
                searchArray = i;
                searchN = 1; % search array size
                currInd = 1; % current index in stack
                pathStack = []; % store the search path index in searchArray
                pathN =  0; % path stack size
                while searchN > 0
                    currentNode = searchArray(currInd); % current node in DFS
                    nextNodes = G(G(:,1) == currentNode,2); % all nodes having an edge from current
                    nextNodes = nextNodes(~ex(nextNodes)); % among which find the nodes not explored
                    
                    if ~isempty(nextNodes) % if any node is found
                        % mark them as explored and set their leaders as source
                        ex(nextNodes) = 1;
                        leader(nextNodes) = s;
                        
                        % push the current index to path
                        pathN = pathN + 1;
                        pathContent = numel(pathStack);
                        if pathN > pathContent
                            copyStack = zeros(1,2*pathContent);
                            copyStack(1:pathN-1) = pathStack;
                            pathStack = copyStack;
                        end
                        pathStack(pathN) = currInd;
                        
                        % current search index go to the first of nextNodes
                        currInd = searchN + 1;
                        
                        % push nextNodes onto searchArray
                        for ne = 1:numel(nextNodes)
                            searchN = searchN + 1;
                            searchContent = numel(searchArray);
                            if searchN > searchContent
                                copyStack = zeros(1,2*searchContent);
                                copyStack(1:searchN-1) = searchArray;
                                searchArray = copyStack;
                            end
                            searchArray(searchN) = nextNodes(ne);
                        end
                        
                    else % if not found
                        % store finishing time to current node
                        t = t + 1;
                        f(currentNode)  = t;
                        
                        % search the next node in nextNodes of the previous node
                        if currInd < searchN
                            searchArray(currInd:(searchN-1)) = searchArray((currInd+1):searchN); % remove currentNode from searchArray
                        else % pop pathStack to go to previous node
                            if pathN > 0
                                currInd = pathStack(pathN);
                            else
                                currInd = 0;
                            end
                            pathN = pathN - 1;
                            pathContent = numel(pathStack);
                            if pathN < pathContent / 4 && pathContent >= 4
                                pathStack = pathStack(1:(pathContent/2));
                            end
                        end
                        
                        % reduce searchArray
                        searchN = searchN - 1;
                        searchContent = numel(searchArray);
                        if searchN < searchContent / 4 && searchContent >= 4
                            searchArray = searchArray(1:(searchContent/2));
                        end
                        
                    end
                    
                end
                
            end
        end
    end

end