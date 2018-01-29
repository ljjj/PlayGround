function [f, x] = wepdf(y, binEdges, w)
% weighted empirical PDF at specified bin edges

    N = numel(y);
    if ~exist('w', 'var')
        w = ones(N,1);
    end
    
    [y, id] = sort(y);
    w = w(id);
    
    binNumber = numel(binEdges) - 1;
    f = zeros(1,binNumber);
    j = 1;
    for i = 1:binNumber % hist counts
        if y(j) > binEdges(i+1) % no data in this bin
            continue
        end
        k = 0;
        while y(j + k) < binEdges(i+1) % closed-open interval [)
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
    binCtrs = (binEdges(1:end-1) + binEdges(2:end))/2;
    f = f/trapz(binCtrs, f); % normalize pdf
    x = (binEdges(1:end-1) + binEdges(2:end))/2;
end