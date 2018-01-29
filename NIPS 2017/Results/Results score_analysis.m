%% read and preprocess data
% load raw pairwise scores
hit_target_1 = csvread('1_hit_target_class_matrix.csv',1,1);
hit_target_2 = csvread('2_hit_target_class_matrix.csv',1,1);
err_1 = csvread('1_error_matrix.csv',1,1);
err_2 = csvread('2_error_matrix.csv',1,1);
acc_2 = csvread('2_accuracy_matrix.csv',1,1);
N_images = 1000;

% load team ids
[n_dfs_2, n_atks_2] = size(err_2);
n_ntg_atks_2 = 62;
dfs_ids_2 = csvread('2_error_matrix.csv',1,0,[1 0 n_dfs_2 0]);
ntg_atks_ids_2 = csvread('2_error_matrix.csv',0,1,[0 1 0 n_ntg_atks_2]);
tg_atks_ids_2 = csvread('2_error_matrix.csv',0,n_ntg_atks_2+1,[0 n_ntg_atks_2+1 0 n_atks_2]);

% load leaderboard results
T_dfs_2 = readtable('2_defense_results.csv');
T_ntg_atks_2 = readtable('2_non_targeted_attack_results.csv');
T_tg_atks_2 = readtable('2_targeted_attack_results.csv');
%% definitions
ftsz = 10; % font size to show team names
select_team_name = false;
if select_team_name
    % select team names highlights
    dfs_1_rows = [49 50 51];
    dfs_1_names = {'bl adv inc v3'...
        'bl ens adv inc res v2'...
        'bl inc v3'};
    
    ntg_atks_1_cols = [28 55 56 57];
    ntg_atks_1_names = {'yaozhao' 'bl fgsm' 'bl noop'...
        'bl randnoise'};
    
    tg_atks_1_cols = [86 87 88];
    tg_atks_1_names = {'bl itertc 10' 'bl itertc 20'...
        'bl steptc'};
    
    
    dfs_2_rows = [2 18 20 23 28 ...
        29 32 36 40 ...
        48 58 59 62 ...
        68 ...
        69 70 71];
    dfs_2_names = {'confirm' 'tonyyy' 'luizgh' 'yaozhao' 'toshik'...
        'rwightman' 'scottclowe' 'mzweilin' 'sunnyraj10411'...
        'rafaelmm' 'ndas30' 'cihangxie' 'hughbzhang' ...
        'aditirag' ...
        'bl adv inc v3' 'bl ens adv inc res v2' 'bl inc v3'};
    
    ntg_atks_2_cols = [2 3 5 10 ...
        12 17 24 ...
        27 28 31 37 39 ...
        38 41 43 47 57 ...
        60 61 62];
    ntg_atks_2_names = {'confirm' 'carpedm20' 'luizgh' 'erwoka' ...
        'taesikna' 'rafaelmm' 'yaozhao' ...
        'voilin' 'rickychou' 'oavdeev' 'mike1808' 'dtsipras' ...
        'sangxia' '-rwightman' 'takeshimiyauchi' '-yinpeng' 'toshik' ...
        'bl fgsm' 'bl noop' 'bl randnoise'};
    
    tg_atks_2_cols = [64 66 68 69 ...
        73 75 76 78 ...
        79 82 85 87 ...
        90 91 94 ...
        95 102 104 ...
        107 108 109];
    tg_atks_2_names = {'confirm' 'erkowa' 'tarobxl' 'shujian' ...
        'luizgh' 'anlthms' 'yaozhao' 'voilin' ...
        'rickychou' '-scottclowe' 'sibosong' '-yinpeng' ...
        '-rwightman' 'mike1808' 'dtsipras' ...
        'stevenxtw' 'shijing' 'lujing' ...
        'bl itertc 10' 'bl itertc 20' 'bl steptc'};
else
    total_2 = acc_2+err_2;
    
    % fill all team names
    dfs_2_rows = nan(size(dfs_ids_2));
    [~, i_dfs_ids_2] = sort(dfs_ids_2);
    [~, i_T_dfs_ids_2] = sort(T_dfs_2.KaggleTeamId);
    dfs_2_rows(i_T_dfs_ids_2) = i_dfs_ids_2;
    dfs_2_names = T_dfs_2.TeamName;
    dfs_2_MET = T_dfs_2.MedianEvalTime;
    row_incomplete = find(max(total_2,[],2) < N_images);
    for j = row_incomplete
        indj = find(dfs_2_rows==j);
        dfs_2_names(indj) = {[num2str(max(total_2(j,:))) ' ' dfs_2_names{indj}]};
    end
    for j = 1:numel(dfs_2_names)
        if T_dfs_2.NormalizedScore(j) < 0
            dfs_2_names(j) = {['---' dfs_2_names{j}]};
        end
    end
    
    ntg_atks_2_cols = nan(size(ntg_atks_ids_2));
    [~, i_ntg_atks_ids_2] = sort(ntg_atks_ids_2);
    [~, i_T_ntg_atks_ids_2] = sort(T_ntg_atks_2.KaggleTeamId);
    ntg_atks_2_cols(i_T_ntg_atks_ids_2) = i_ntg_atks_ids_2;
    ntg_atks_2_names = T_ntg_atks_2.TeamName;
    ntg_atks_2_MET = T_ntg_atks_2.MedianEvalTime;
    col_incomplete = find(max(total_2) < N_images);
    for j = col_incomplete
        indj = find(ntg_atks_2_cols==j);
        if indj <= n_ntg_atks_2
            ntg_atks_2_names(indj) = {[num2str(max(total_2(:,j))) ' ' ntg_atks_2_names{indj}]};
        end
    end
    for j = 1:numel(ntg_atks_2_names)
        if T_ntg_atks_2.NormalizedScore(j) < 0
            ntg_atks_2_names(j) = {['---' ntg_atks_2_names{j}]};
        end
    end
    
    tg_atks_2_cols = nan(size(tg_atks_ids_2));
    [~, i_tg_atks_ids_2] = sort(tg_atks_ids_2);
    [~, i_T_tg_atks_ids_2] = sort(T_tg_atks_2.KaggleTeamId);
    tg_atks_2_cols(i_T_tg_atks_ids_2) = n_ntg_atks_2 + i_tg_atks_ids_2;
    tg_atks_2_names = T_tg_atks_2.TeamName;
    tg_atks_2_MET = T_tg_atks_2.MedianEvalTime;
    for j = col_incomplete
        indj = find(tg_atks_2_cols==j);
        if indj > n_ntg_atks_2
            tg_atks_2_names(indj) = {[num2str(max(total_2(:,j))) ' ' tg_atks_2_names{indj}]};
        end
    end
    for j = 1:numel(tg_atks_2_names)
        if T_tg_atks_2.NormalizedScore(j) < 0
            tg_atks_2_names(j) = {['---' tg_atks_2_names{j}]};
        end
    end
end
%% plot h1
figure
b = bar3(hit_target_1);
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
hold on
for i = 1:numel(dfs_1_rows)
    text(-10, dfs_1_rows(i), 1000,...
        dfs_1_names{i}, 'fontsize', ftsz)
end
for i = 1:numel(ntg_atks_1_cols)
    text(ntg_atks_1_cols(i), 0, 1000,...
        ntg_atks_1_names{i}, 'rotation', 90, 'fontsize', ftsz)
end
for i = 1:numel(tg_atks_1_cols)
    text(tg_atks_1_cols(i), 0, 1000,...
        tg_atks_1_names{i}, 'rotation', 90, 'fontsize', ftsz)
end
title('hit target class matrix 1')
axis tight
view(2)
%% plot h2
figure
b = bar3(hit_target_2);
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
hold on
plot([62.5 62.5], [0.5 71.5], 'r', 'lineWidth', 3)
for i = 1:numel(dfs_2_rows)
    text(-10, dfs_2_rows(i), 1000,...
        dfs_2_names{i}, 'fontsize', ftsz)
end
for i = 1:numel(ntg_atks_2_cols)
    text(ntg_atks_2_cols(i), 0, 1000,...
        ntg_atks_2_names{i}, 'rotation', 90, 'fontsize', ftsz)
end
for i = 1:numel(tg_atks_2_cols)
    text(tg_atks_2_cols(i), 0, 1000,...
        tg_atks_2_names{i}, 'rotation', 90, 'fontsize', ftsz)
end
title('hit target class matrix 2')
axis tight
view(2)
%% plot e1
figure
b = bar3(err_1);
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
hold on
for i = 1:numel(dfs_1_rows)
    text(-10, dfs_1_rows(i), 1000,...
        dfs_1_names{i}, 'fontsize', ftsz)
end
for i = 1:numel(ntg_atks_1_cols)
    text(ntg_atks_1_cols(i), 0, 1000,...
        ntg_atks_1_names{i}, 'rotation', 90, 'fontsize', ftsz)
end
for i = 1:numel(tg_atks_1_cols)
    text(tg_atks_1_cols(i), 0, 1000,...
        tg_atks_1_names{i}, 'rotation', 90, 'fontsize', ftsz)
end
title('error matrix 1')
axis tight
view(2)
%% plot e2
figure
b = bar3(err_2);
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
hold on
plot([62.5 62.5], [0.5 71.5], 'r', 'lineWidth', 3)
for i = 1:numel(dfs_2_rows)
    text(-10, dfs_2_rows(i), 1000,...
        dfs_2_names{i}, 'fontsize', ftsz)
end
for i = 1:numel(ntg_atks_2_cols)
    text(ntg_atks_2_cols(i), 0, 1000,...
        ntg_atks_2_names{i}, 'rotation', 90, 'fontsize', ftsz)
end
for i = 1:numel(tg_atks_2_cols)
    text(tg_atks_2_cols(i), 0, 1000,...
        tg_atks_2_names{i}, 'rotation', 90, 'fontsize', ftsz)
end
title('error class matrix 2')
axis tight
view(2)
%% sort err_2
% err_2 = hit_target_2;
% sort rows first for both non-target and target attacks agains the same defense
custom_sort = false;
if custom_sort
    [~,ie2r] = sort(err_2(:,9));
else
    [~,i_dfs_score_2] = sort(abs(T_dfs_2.NormalizedScore));
    ie2r = nan(size(dfs_ids_2));
    [~, i_dfs_ids_2] = sort(dfs_ids_2);
    [~, i_T_dfs_ids_score_2] = sort(T_dfs_2.KaggleTeamId(i_dfs_score_2));
    ie2r(i_T_dfs_ids_score_2) = i_dfs_ids_2;
end
err_2_sorted = err_2(ie2r,:);
hit_target_2_sorted = hit_target_2(ie2r,:);

% sort columns for non-target and target attacks separately
err_2_ntg = err_2_sorted(:,1:n_ntg_atks_2);
% err_2_tg = err_2_sorted(:,n_ntg_atks_2+1:end);
err_2_tg = hit_target_2_sorted(:,n_ntg_atks_2+1:end);
if custom_sort
    [~,ie2n] = sort(err_2_ntg(70,:));
    [~,ie2t] = sort(err_2_tg(70,:));
    % ie2n = 1:size(err_2_ntg,2);
    % ie2t = 1:size(err_2_tg,2);
else
    [~,i_ntg_atks_score_2] = sort(abs(T_ntg_atks_2.NormalizedScore));
    ie2n = nan(size(ntg_atks_ids_2));
    [~, i_ntg_atks_ids_2] = sort(ntg_atks_ids_2);
    [~, i_T_ntg_atks_ids_score_2] = sort(T_ntg_atks_2.KaggleTeamId(i_ntg_atks_score_2));
    ie2n(i_T_ntg_atks_ids_score_2) = i_ntg_atks_ids_2;
    
    [~,i_tg_atks_score_2] = sort(abs(T_tg_atks_2.NormalizedScore));
    ie2t = nan(size(tg_atks_ids_2));
    [~, i_tg_atks_ids_2] = sort(tg_atks_ids_2);
    [~, i_T_tg_atks_ids_score_2] = sort(T_tg_atks_2.KaggleTeamId(i_tg_atks_score_2));
    ie2t(i_T_tg_atks_ids_score_2) = i_tg_atks_ids_2;
end
err_2_sorted = [err_2_ntg(:,ie2n) err_2_tg(:,ie2t)];

% plot sorted e2
figure('visible','off')
b = bar3(err_2_sorted);
c = colorbar;
title(c,'error matrix 2')
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
hold on
plot(n_ntg_atks_2 + [0.5 0.5], [0.5 n_dfs_2+0.5], 'r', 'lineWidth', 3)
plot(n_atks_2 + [0.5 0.5], [0.5 n_dfs_2+1.5], 'r', 'lineWidth', 3)
plot([0.5 n_atks_2+1.5], n_dfs_2 + [0.5 0.5], 'r', 'lineWidth', 3)

dfs_2_MET_sorted = zeros(size(dfs_2_MET));
for i = 1:numel(dfs_2_rows)
    new_id = find(ie2r==dfs_2_rows(i));
    text(-10, new_id, 1000, dfs_2_names{i}, 'fontsize', ftsz)
    dfs_2_MET_sorted(new_id) = dfs_2_MET(i);
end
text(-10, n_dfs_2+1, 1000, '2x ATK MED TIME', 'fontsize', ftsz)

ntg_atks_2_MET_sorted = zeros(size(ntg_atks_2_MET));
for i = 1:numel(ntg_atks_2_cols)
    new_id = find(ie2n==ntg_atks_2_cols(i));
    text(new_id, 0, 1000, ntg_atks_2_names{i}, 'rotation', 90, 'fontsize', ftsz)
    ntg_atks_2_MET_sorted(new_id) = ntg_atks_2_MET(i);
end

tg_atks_2_MET_sorted = zeros(size(tg_atks_2_MET));
for i = 1:numel(tg_atks_2_cols)
    new_id = find(ie2t==tg_atks_2_cols(i)-n_ntg_atks_2);
    text(n_ntg_atks_2+new_id, 0, 1000, tg_atks_2_names{i}, 'rotation', 90, 'fontsize', ftsz)
    tg_atks_2_MET_sorted(new_id) = tg_atks_2_MET(i);
end
text(n_atks_2+1, 0, 1000, '2x DFS MED TIME', 'rotation', 90, 'fontsize', ftsz)

b = bar3([err_2_sorted dfs_2_MET_sorted*2; ...
    ntg_atks_2_MET_sorted'*2 tg_atks_2_MET_sorted'*2 0]);
c = colorbar;
caxis([0 1000])
title(c,'error matrix 2')
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

% title('error matrix 2')
axis tight
view(2)

% save plot
set(gcf, 'PaperPosition', [0 0 1 0.8] * 20);
saveName = 'Non-target Error v Target Hit Matrix 2';
% saveName = 'error matrix 2';
% saveName = 'hit target matrix 2';
disp(['Saving file ', saveName,' ...'])
print(gcf, [saveName,'.jpg'], '-djpeg', '-r300', '-opengl')
close(gcf)