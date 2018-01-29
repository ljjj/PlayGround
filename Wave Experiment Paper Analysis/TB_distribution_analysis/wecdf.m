function [f, x] = wecdf(y, bin_edges, w)
% weighted empirical CDF at specified bin edges

    N = numel(y);
    if ~exist('w', 'var')
        w = ones(N,1);
    end
    
    [y, id] = sort(y);
    w = w(id);
    
    bin_number = numel(bin_edges) - 1;
    f = zeros(1,bin_number);
    j = 1;
    for i = 1:bin_number % hist counts
        if y(j) > bin_edges(i+1) % no data in this bin
            continue
        end
        k = 0;
        while y(j + k) < bin_edges(i+1) % closed-open interval [)
            k = k + 1;
            if j + k > N 
                break
            end
        end
        f(i) = sum(w(j:(j+k-1)));
        j = j + k;
        if j > N
            break
        end
    end
    f(end) = f(end) + numel(y(j:end));
    f = cumsum(f); % cumulative counts
    f = f/f(end); % cumulative probability
    x = (bin_edges(1:end-1) + bin_edges(2:end))/2;
end