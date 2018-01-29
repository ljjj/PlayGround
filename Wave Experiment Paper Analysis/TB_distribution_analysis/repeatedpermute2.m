function [x1, y1] = repeatedpermute2(x,y)
% permute labels between 2 samples with contiguous repeats always together

% standardize input
x = x(:);
y = y(:);

% reduce the repeats
[ux,rx] = repeatreduce(x);
[uy,ry] = repeatreduce(y);

% shuffle values
Nx = numel(ux);
Ny = numel(uy);
uz = [ux; uy];
rz = [rx; ry];
shuffle_ind = randperm(Nx+Ny);
uz = uz(shuffle_ind);
rz = rz(shuffle_ind);
ux1 = uz(1:Nx);
uy1 = uz(Nx+1:Nx+Ny);
rx1 = rz(1:Nx);
ry1 = rz(Nx+1:Nx+Ny);

% expand the repeats
x1 = repeatexpand(ux1,rx1);
y1 = repeatexpand(uy1,ry1);
end