% This file calculates wave positions using Gaussmian mixture fits.
%% load data
load('..\..\..\fig2_rawdata.mat')
load('..\Gaussian_mixture_fits.mat')
load('..\wave_speeds_and_cell_counts.mat')
Nexp = numel(data2);
count_thresh = cell(1, Nexp);
%% show plots
wave_end_pos = cell(1, Nexp);
wave_end_pos_b = cell(1, Nexp);
wave_end_pos_f = cell(1, Nexp);
wave_end_time = cell(1, Nexp);
wave_end_time_b = cell(1, Nexp);
wave_end_time_f = cell(1, Nexp);
for i = 1:Nexp % different experiments
    % initialization
    if isempty(data2(i).CellDensity)
        continue
    end
    t_str = ['Asp = ',num2str(data2(i).AspConc),' \muM'];
    if isempty(data2(i).TBdist)
        t_str = [t_str, ' no P(TB)'];
    else
        t_str = [t_str, ' has P(TB)'];
    end
    cp = data2(i).CellDensity.CellPos;
    if isempty(cp(2).pos)
        sweeps = 1:2:numel(cp);
    else
        sweeps = 1:numel(cp);
    end
    N_s = numel(sweeps);
    
    % find approximate number of cells in the wave
    cm = cp(5).counts_smooth; % later wave for accuracy but not too late to keep more cells
    dsd = diff(sign( diff(cm) ));
    den_min_ind = find(dsd == 2) + 1; % find where first derivative = 0 and second derivative > 0
    if isempty(den_min_ind)
        den_min_ind = 1;
    end
    den_min_ind = den_min_ind(1);
    cc = cumsum(cp(5).counts, 'reverse');
    total_count = cc(den_min_ind);
    % used to determine the maximal threshold to use
    count_max = round(total_count*1.6);
    dc = round(count_max/300);
    count_thresh{i} = dc:dc:count_max; % use cell count threshold to determine wave front
    N_t = numel(count_thresh{i});
    th_pos = nan(N_s, N_t);
    th_time = nan(N_s, N_t);
    
    for j = 1:N_s % different sweeps
        cum_counts = cumsum(cp(sweeps(j)).counts, 'reverse');
        for k = 1:N_t % different wave front threshold choices
            ind = find(cum_counts > count_thresh{i}(k));
            if isempty(ind)
                l = 1;
            else
                l = ind(end); % index where cell count exceeds threshold
            end
            th_pos(j,k) = cp(sweeps(j)).pos(l);
            th_time(j,k) = cp(sweeps(j)).time(l);
        end
    end

    % index to define where the wave is
    wave_end = speeds_ls{i} < mus(i,2);
    wave_end_l = speeds_ls{i} < mus(i,2) - sigmas(i,2);
    wave_end_h = speeds_ls{i} < mus(i,2) + sigmas(i,2);
    [~, imax] = max(speeds_ls{i});
    
    % not in the front of the wave
    wave_end(1:imax) = 0;
    wave_end_l(1:imax) = 0;
    wave_end_h(1:imax) = 0;
    
    for j = 1:numel(speeds_ls{i})
        if wave_end(j)
            wave_end_pos{i} = th_pos(:,j);
            wave_end_time{i} = th_time(:,j);
            break
        end
    end
    
    for j = 1:numel(speeds_ls{i})
        if wave_end_l(j)
            wave_end_pos_b{i} = th_pos(:,j);
            wave_end_time_b{i} = th_time(:,j);
            break
        end
    end
    
    for j = 1:numel(speeds_ls{i})
        if wave_end_h(j)
            wave_end_pos_f{i} = th_pos(:,j);
            wave_end_time_f{i} = th_time(:,j);
            break
        end
    end
end
%% save final results
save('wave_ends.mat','wave_end_pos','wave_end_pos_b','wave_end_pos_f','wave_end_time','wave_end_time_b','wave_end_time_f')