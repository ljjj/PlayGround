%% load
dropbox_dir = 'C:\Users\jl2345\Dropbox (emonetlab)\';
% dropbox_dir = 'D:\Dropbox (emonetlab)\';
raw_data_file = [dropbox_dir, 'paper_wave_diversity\manuscript\Figures\final figure May2016\fig2\fig2_material_May2017\fig2_rawdata.mat'];
load(raw_data_file)
%% calculate quantiles
tb_q = nan(3,3,2); % tb quantiles
for e = 1:3
    tb_q(e,1,:) = quantile(data2(3*e-2).TB, [0.7 0.8]);
    tb_q(e,2,:) = quantile(data2(3*e-1).TB, [0.7 0.8]);
    tb_q(e,3,:) = quantile(data2(3*e).TB, [0.7 0.8]);
end
mean_tb_q = squeeze(mean(tb_q, 2));
std_tb_q = squeeze(std(tb_q, 0, 2));
save('tb_quantiles', 'mean_tb_q', 'std_tb_q', 'tb_q')