function y = repeatedsample(x)
% randomly sample with contiguous repeats always together, resulting in the
% same number of sample blocks

% standardize input
x = x(:);

% reduce the repeats
[ux,rx] = repeatreduce(x);

% sample data
Nx = numel(ux);
[uy, idx] = datasample(ux, Nx, 'Weights', rx);
ry = rx(idx);

% expand the repeats
y = repeatexpand(uy,ry);
end