function v = repeatexpand(u,r)
% u     array of unique entries
% r     number of repeats for each entry
% v     array with entries repeated
    u = u(:);
    r = r(:);
    Nr = numel(u);
    Nv = sum(r);
    v = nan(Nv, 1);
    cr = [0; cumsum(r)];
    for i = 1:Nr
        v(cr(i)+1:cr(i+1)) = u(i);
    end
end