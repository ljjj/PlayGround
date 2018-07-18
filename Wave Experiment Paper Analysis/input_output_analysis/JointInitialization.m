function Joint = JointInitialization(Series1, Series2, dt, xc_max_lag, mi_max_lag, lf_max_lag)
% This function initializes Joint struct for certain calculations between
% two Series structs.
% Joint : struct with fields:
%         xc_max_lag : max lag of cross-correlation between df/dt and dx/dt
%         xc         : cross-correlation values
%         xc_lags    : cross-correlation lags
%                      xc_lags < 0 means Series1 causes Series2
%         mi_max_lag : max lag to calculate mutual information
%         mi         : mutual information values
%         mi_lags    : mutual information lags
%                      mi_lags < 0 means Series1 causes Series2
%         lf_max_lag : max lag to calculate linear filter
%         lf         : linear filter values
%         lf_lags    : linear filter lags
%                      lf_lags > 0 means Series1 causes Series2
%         lf_r2      : goodness of fit of the linear filter

Joint = [];
Joint.name1 = Series1.name;
Joint.name2 = Series2.name;
Joint.xc_max_lag = xc_max_lag;
Joint.xc = [];
Joint.xc_lags = (-xc_max_lag:xc_max_lag)*dt;
Joint.mi_max_lag = mi_max_lag;
Joint.mi = [];
Joint.mi_lags = (-mi_max_lag:mi_max_lag)*dt;
Joint.lf_max_lag = lf_max_lag;
Joint.lf = [];
Joint.lf_lags = (0:lf_max_lag)*dt;
Joint.lf_r2 = [];
end