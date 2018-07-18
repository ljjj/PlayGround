function Joint = JointAppend(Joint, Series1, Series2)
% This function appends Joint calculations of two newly appended Series
% structs to Joint struct.
% Joint     : struct defined in JointInitialization.m
% Series1,2 : structs defined in SeriesInitialization.m

% check input Series consistency
if ~strcmp(Joint.name1, Series1.name)
    error('Series 1 mismatch')
end
if ~strcmp(Joint.name2, Series2.name)
    error('Series 2 mismatch')
end

% append joint calculations
ts1 = Series1.timeseries(end,:);
ts2 = Series2.timeseries(end,:);

% cross correlation
xc = xcorr(ts1, ts2, Joint.xc_max_lag, 'coeff');
Joint.xc = [Joint.xc; xc(:)'];

% mutual information
mi = nan(1, 2*Joint.mi_max_lag+1);
for l = -Joint.mi_max_lag:Joint.mi_max_lag
    mi(l+Joint.mi_max_lag+1) = SeriesMutualInformation(Series1, Series2, l);
end
Joint.mi = [Joint.mi; mi];

% linear filter
x = [ts1(:); zeros(Joint.lf_max_lag,1)];
X = x;
for j = 1:Joint.lf_max_lag
    X = [X lagmatrix(x,j)];
end
X(isnan(X)) = 0;
y = [ts2(:); zeros(Joint.lf_max_lag,1)];
lf = pinv(X'*X)*X'*y;
yp = fftfilt(lf, ts1);
lf_r2 = 1-sum((yp(:) - ts2(:)).^2)/sum(ts2.^2);
Joint.lf = [Joint.lf; lf'];
Joint.lf_r2 = [Joint.lf_r2; lf_r2];
end