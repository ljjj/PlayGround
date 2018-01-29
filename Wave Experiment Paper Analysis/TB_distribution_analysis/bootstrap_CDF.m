%% load
data_choice = 1;
bs_method = 2; % 1 = frame-wise, 2 = track-wise
dropbox_dir = 'C:\Users\jl2345\Dropbox (emonetlab)\';
% dropbox_dir = 'D:\Dropbox (emonetlab)\';
raw_data_file = [dropbox_dir, 'paper_wave_diversity\manuscript\Figures\codes for figs\raw data\data',num2str(data_choice),'.mat'];
load(raw_data_file)
if data_choice == 1
    data1(11) = data1(9);
    data1(9) = data1(10);
    data1(10) = data1(11);
    data1(11) = [];
    data = data1;
    clear data1
    N_group = 3; % [wave1, wave2, ori]
elseif data_choice == 2
    data = data2;
    clear data2
    N_group = 4; % [50, 100, 200, 0]
end
N_replicates = 3; % number of replicates in each condition
N_exps = N_group*N_replicates; % total number of experiments
%% bootstrap for CDF error
bin_number = 50;
bin_edges = 0:1/bin_number:1;
bin_ctrs = (bin_edges(1:end-1) + bin_edges(2:end))/2;
bs_cdf = cell(1,N_exps);
N_bootstrap = 10000;
for e = 1:N_exps
    disp(['Boostrapping experiment ', num2str(e)])
    if bs_method == 1
        bs_cdf{e} = bootstrp(N_bootstrap, @(x) wecdf(x, bin_edges), data(e).TB); % frame-wise bootstrapping
    elseif bs_method == 2
        % track-wise bootstrapping
        bs_cdf{e} = zeros(N_bootstrap,bin_number);
        x = data(e).TB(:);
        [ux,rx] = repeatreduce(x);
        Nx = numel(ux);
        for b = 1:N_bootstrap
            % sample data
            [uy, idx] = datasample(ux, Nx, 'Weights', rx);
            ry = rx(idx);
            bs_cdf{e}(b,:) = wecdf(uy, bin_edges, ry);
        end
    end
end
%% calculate CDF error from bootstrap
bs_cdf_mean = nan(N_exps, bin_number);
bs_cdf_std = nan(N_exps, bin_number);
for e = 1:N_exps
    bs_cdf_mean(e,:) = mean(bs_cdf{e});
    bs_cdf_std(e,:) = std(bs_cdf{e});
end

bs_cdf_std_max = nan(1,N_exps);
for e = 1:N_exps
    bs_cdf_std_max(e) = max(bs_cdf_std(e,:));
%     disp(['Max std in CDF of experiment ', num2str(e), ' is ', num2str(bs_cdf_std_max(e))])
end

bs_cdf_abs_diff_max = nan(N_exps, N_exps);
for e1 = 1:N_exps
    for e2 = 1:N_exps
        bs_cdf_abs_diff_max(e1, e2) = max(abs(bs_cdf_mean(e1,:) - bs_cdf_mean(e2,:)));
%         disp(['Max abs diff in CDF of experiment ', num2str(e1), ' and experiment ',...
%             num2str(e2), ' is ', num2str(bs_cdf_abs_diff_max(e1, e2))])
    end
end
%% save boostrap results
if bs_method == 1
    save(['bootstrapped_CDF_frame-wise_data',num2str(data_choice),'.mat'], 'bin_number', 'bin_edges', 'bin_ctrs',...
        'bs_cdf', 'bs_cdf_mean', 'bs_cdf_std', 'bs_cdf_std_max', 'bs_cdf_abs_diff_max')
elseif bs_method == 2
    save(['bootstrapped_CDF_track-wise_data',num2str(data_choice),'.mat'], 'bin_number', 'bin_edges', 'bin_ctrs',...
        'bs_cdf', 'bs_cdf_mean', 'bs_cdf_std', 'bs_cdf_std_max', 'bs_cdf_abs_diff_max')
end
%% load results
if bs_method == 1
    load(['bootstrapped_CDF_frame-wise_data',num2str(data_choice),'.mat'])
elseif bs_method == 2
    load(['bootstrapped_CDF_track-wise_data',num2str(data_choice),'.mat'])
end
%% show CDF differences colorbar
figure
b = bar3(bs_cdf_abs_diff_max);
c = colorbar;
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
xlabel('Experiment #')
ylabel('Experiment #')
caxis([0 0.4])
zlim([0 0.45])
title('sup|CDF(i)-CDF(j)|')