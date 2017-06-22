function result = permutationtest2(x, y, k)
% x and y are two samples from two distributions
% we permute labels to test whether they are two different distributions

% find the two sample KS statistics
[~,~,kstat0] = kstest2(x,y);

% permute sample labels and find KS statistics
kstats = zeros(1,k);
for i = 1:k
    disp(['permuting ', num2str(i), ' out of ', num2str(k)])
    [x1, y1] = repeatedpermute2(x,y);
    [~,~,kstati] = kstest2(x1,y1);
    kstats(i) = kstati;
end

% find p-value from permuted KS statistics
kstats = sort(kstats);
j = 1;
while kstats(j) < kstat0
    j = j + 1;
    if j > k
        break
    end
end
j = j - 1;
p = 1 - j/k;

result = [p kstat0 kstats];
end