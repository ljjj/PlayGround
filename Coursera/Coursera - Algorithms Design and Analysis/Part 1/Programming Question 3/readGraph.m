function graph = readGraph(file)
    fid = fopen(file);
    
    edgeID = 0;
    graph = [0 0];
    
    % read line by line
    line = fgets(fid);
    while ischar(line)
        edges = str2num(line);
        % the first element in the line is one vertex while the remaining
        % elements are the vertices of all the edges connected to it
        for j = 2:numel(edges)
            % only add edge if not already added
            if ~any(graph(graph(:,1) == edges(j),2) == edges(1))
                edgeID = edgeID + 1;
                graph(edgeID,:) = [edges(1) edges(j)];
            end
        end
        line = fgets(fid);
    end
    
    fclose(fid);
end