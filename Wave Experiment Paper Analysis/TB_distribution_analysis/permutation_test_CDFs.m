%% load and initialize
data_choice = 1;
dropbox_dir = 'C:\Users\jl2345\Dropbox (emonetlab)\';
% dropbox_dir = 'D:\Dropbox (emonetlab)\';
raw_data_file = [dropbox_dir, 'paper_wave_diversity\manuscript\Figures\codes for figs\raw data\data',num2str(data_choice),'.mat'];
load(raw_data_file)

if data_choice == 1
    % fix data order
    data1(11) = data1(9);
    data1(9) = data1(10);
    data1(10) = data1(11);
    data1(11) = [];
    data = data1;
    clear data1
    N_group = 3; % number of conditions
    % pairwise: [wave1&wave2, wave2&ori, ori&wave1]
    group_str = {'wave1', 'wave2', 'ori'}; % name of each group
    colors = {'r', [1 0.843137254901961 0], 'k'}; % color of each group
elseif data_choice == 2
    data = data2;
    clear data2
    N_group = 4; % number of conditions
    % pairwise: [200&0, 100&0, 100&200, 50&0, 50&200, 50&100]
    group_str = {'50 \muM', '100 \muM', '200 \muM', 'WT'}; % name of each group
    colors = {'r', [1 0.843137254901961 0], 'b', 'k'}; % color of each group
end

N_replicates = 3; % number of replicates in each condition
N_exps = N_group*N_replicates; % total number of experiments

in_group_inds = combnk(1:N_replicates,2); % indices for within group pairwise comparisons
in_group_inds_1 = in_group_inds(:,1);
in_group_inds_2 = in_group_inds(:,2);

out_group_inds = combnk(1:N_group,2); % indices for pairwise group comparisons
out_group_inds_1 = out_group_inds(:,1);
out_group_inds_2 = out_group_inds(:,2);

repicate_label_permutations = unique(perms([ones(N_replicates,1) 2*ones(N_replicates,1)]), 'rows'); % indices for permutations of replicate labels in two groups
N_replicate_permutations = nchoosek(N_replicates*2, N_replicates)/2; % number of permutations of the replicates in two groups
repicate_label_permutations = repicate_label_permutations(1:N_replicate_permutations,:);

group_choices = combnk(1:(N_group*N_replicates), 2*N_replicates); % indices for choices of pairwise group comparisons
N_all_permutations = nchoosek(N_group*N_replicates,2*N_replicates); % number of choices of two groups among all replicates
%% show CDF plots
bin_number = 50;
bin_edges = 0:1/bin_number:1;
bin_ctrs = (bin_edges(1:end-1) + bin_edges(2:end))/2;
figure
for g = 1:N_group
    ex = ((g-1)*N_replicates+1):(g*N_replicates);
    x = vertcat(data(ex).TB);
    [nt,ct] = wecdf(x,bin_edges);
    subplot(2,1,1)
    hold on
    plot(ct, nt, 'Color', colors{g}, 'LineWidth', 3)
    for e = ex
        [n,c] = wecdf(data(e).TB,bin_edges);
        plot(c, n, 'Color', colors{g}, 'LineStyle', '--')
    end
    xlabel('TB')
    ylabel('CDF')
    subplot(2,1,2)
    hold on
    plot(ct, diff([0 nt]), 'Color', colors{g}, 'LineWidth', 3)
    for e = ex
        [n,c] = wecdf(data(e).TB,bin_edges);
        plot(c, diff([0 n]), 'Color', colors{g}, 'LineStyle', '--')
    end
    xlabel('TB')
    ylabel('PDF')
end
%% calculate p value pairwise between different distributions using permutation test
N_perms = 1000; % number of MC permutations

tb_dist_p_in_group = nan(N_group, N_replicates, N_perms+2); % groups x [1&2, 2&3, 3&1]
for group = 1:N_group
    parfor in_group_ind = 1:N_replicates
        e1 = in_group_inds_1(in_group_ind);
        e2 = in_group_inds_2(in_group_ind);
        ex1 = (group-1)*N_replicates + e1;
        ex2 = (group-1)*N_replicates + e2;
        disp(['Permutation testing experiment ', num2str(ex1), ' and experiment ', num2str(ex2)])
        tb_dist_p_in_group(group,in_group_ind,:) = permutationtest2(data(ex1).TB, data(ex2).TB, N_perms);
    end
end

tb_dist_p_out_group = nan(N_group*(N_group-1)/2, N_perms+2);
parfor out_group_ind = 1:(N_group*(N_group-1)/2)
    g1 = out_group_inds_1(out_group_ind);
    g2 = out_group_inds_2(out_group_ind);
    ex1 = ((g1-1)*N_replicates+1):(g1*N_replicates);
    ex2 = ((g2-1)*N_replicates+1):(g2*N_replicates);
    disp(['Permutation testing pooled experiments ', num2str(ex1), ' and pooled experiments ', num2str(ex2)])
    x = vertcat(data(ex1).TB);
    y = vertcat(data(ex2).TB);
    tb_dist_p_out_group(out_group_ind,:) = permutationtest2(x, y, N_perms);
end
%% show KS statistics
figure
for group = 1:N_group
    for in_group_ind = 1:N_replicates
        subplot(N_group,N_replicates,(group-1)*N_replicates+in_group_ind)
        histogram(tb_dist_p_in_group(group,in_group_ind,3:end), 'normalization', 'cdf')
        hold on
        plot(tb_dist_p_in_group(group,in_group_ind,2)*[1 1],[0 1], 'r', 'LineWidth', 3)
        axis([0 0.25 0 1])
        xlabel('KS statistic')
        ylabel('CDF')
        if group == 1
            title([num2str(in_group_inds_1(in_group_ind)), ' vs. ',...
                num2str(in_group_inds_2(in_group_ind))])
        end
        if in_group_ind == 1
            ylabel(group_str{group})
        end
    end
end

figure
for out_group_ind = 1:(N_group*(N_group-1)/2)
    subplot(N_group*(N_group-1)/2,1,out_group_ind)
    histogram(tb_dist_p_out_group(out_group_ind,3:end), 'normalization', 'cdf')
    hold on
    plot(tb_dist_p_out_group(out_group_ind,2)*[1 1],[0 1], 'r', 'LineWidth', 3)
    axis([0 0.25 0 1])
    xlabel('KS statistic')
    ylabel('CDF')
    title([group_str{out_group_inds_1(out_group_ind)},' vs ',...
        group_str{out_group_inds_2(out_group_ind)}])
end
%% replicate label permutation test
tb_dist_p_out_group_replicate = nan(N_group*(N_group-1)/2,N_replicate_permutations);
for out_group_ind = 1:(N_group*(N_group-1)/2)
    g1 = out_group_inds_1(out_group_ind);
    g2 = out_group_inds_2(out_group_ind);
    ex1 = ((g1-1)*N_replicates+1):(g1*N_replicates);
    ex2 = ((g2-1)*N_replicates+1):(g2*N_replicates);
    disp(['Testing permutation of labels from experiments ', num2str(ex1), ' and experiments ', num2str(ex2)])
    ex = [ex1 ex2];
    parfor rpl = 1:N_replicate_permutations
        exr1 = ex(repicate_label_permutations(rpl,:)==1);
        exr2 = ex(repicate_label_permutations(rpl,:)==2);
        disp(['    Calculating KS statistic of experiments ', num2str(exr1), ' and experiments ', num2str(exr2)])
        x = vertcat(data(exr1).TB);
        y = vertcat(data(exr2).TB);
        [~,~,tb_dist_p_out_group_replicate(out_group_ind,rpl)] = kstest2(x,y);
%         [~,~,tb_dist_p_out_group_replicate(out_group_ind,rlp)] = kstest2(x,y,'Tail','larger');
    end
end
%% show replicate label results
figure
b = bar3(tb_dist_p_out_group_replicate);
c = colorbar;
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
xlabel('Permutations')
ylabel('Group Comparisons')
title('sup|CDF(i)-CDF(j)|')
ax = gca;
for out_group_ind = 1:(N_group*(N_group-1)/2)
    ax.YTickLabels{out_group_ind} = [group_str{out_group_inds_1(out_group_ind)},' vs ',...
        group_str{out_group_inds_2(out_group_ind)}];
end
%% replicate label permutation test across all conditions
tb_dist_p_all_replicate = nan(N_all_permutations,N_replicate_permutations);
parfor gc = 1:N_all_permutations
    ex = group_choices(gc,:);
    disp(['Testing permutation of labels from experiments ', num2str(ex)])
    for rpl = 1:N_replicate_permutations
        exr1 = ex(repicate_label_permutations(rpl,:)==1);
        exr2 = ex(repicate_label_permutations(rpl,:)==2);
        disp(['    Calculating KS statistic of experiments ', num2str(exr1), ' and experiments ', num2str(exr2)])
        x = vertcat(data(exr1).TB);
        y = vertcat(data(exr2).TB);
        [~,~,tb_dist_p_all_replicate(gc,rpl)] = kstest2(x,y);
    end
end
%% show all replicates results
% first determine which ones are the actual group labels
og_ex = nan(N_group*(N_group-1)/2,2*N_replicates);
for out_group_ind = 1:(N_group*(N_group-1)/2)
    g1 = out_group_inds_1(out_group_ind);
    g2 = out_group_inds_2(out_group_ind);
    og_ex(out_group_ind,:) = [((g1-1)*N_replicates+1):(g1*N_replicates) ((g2-1)*N_replicates+1):(g2*N_replicates)];
end
% find them in all permutations
all_perm_ind = nan(N_group*(N_group-1)/2,1);
for gc = 1:N_all_permutations
    for out_group_ind = 1:(N_group*(N_group-1)/2)
        if all(group_choices(gc,:) == og_ex(out_group_ind,:))
            all_perm_ind(out_group_ind) = gc;
        end
    end
end
% plot
figure
histogram(tb_dist_p_all_replicate(:), 'normalization', 'cdf')
hold on
for out_group_ind = 1:(N_group*(N_group-1)/2)
    plot(tb_dist_p_all_replicate(all_perm_ind(out_group_ind),1)*[1 1],...
        [0 0.5], 'Color', colors{out_group_inds_1(out_group_ind)}, 'LineWidth', 3)
    plot(tb_dist_p_all_replicate(all_perm_ind(out_group_ind),1)*[1 1],...
        [0.5 1], 'Color', colors{out_group_inds_2(out_group_ind)}, 'LineWidth', 3)
end
axis([0 0.4 0 1])
xlabel('KS statistics')
ylabel('CDF')
%% save permutation test results
save(['P(TB)_permuation_test_data',num2str(data_choice),'.mat'], 'tb_dist_p_in_group',...
    'tb_dist_p_out_group', 'tb_dist_p_out_group_replicate', 'tb_dist_p_all_replicate', '-v7.3')
%% load
load(['P(TB)_permuation_test_data',num2str(data_choice),'.mat'])