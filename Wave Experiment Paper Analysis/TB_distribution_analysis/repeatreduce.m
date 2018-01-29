function [u,r] = repeatreduce(v)
% v     array with contiguous repeated entries
% u     unique entries
% r     number of repeats
    v = v(:);
    dv = diff(v);
    id = find(dv ~= 0);
    id = [0; id(:); numel(v)];
    Nr = numel(id)-1;
    u = nan(Nr, 1);
    r = nan(Nr, 1);
    for i = 1:Nr
        u(i) = v(id(i+1));
        r(i) = id(i+1) - id(i);
    end
end