traj_length = 120;
%% interactive plot
% load(['stats_',num2str(traj_length),'.mat'])

Series_names = {'dfdt', 'dxdt', 'tb1', 'tb21', 'tb41', 'tb81'};
N_series = numel(Series_names);
N_PC = 3;

all_vars = [];
all_names = {};
for s = 1:N_series
    name = Series_names{s};
    all_vars = [all_vars eval([name, '.pca_results.score(:,1:N_PC)'])];
    for p = 1:N_PC
        all_names = [all_names, ['P(', eval([name, '.name']),') PC',num2str(p)]];
    end
    all_vars = [all_vars eval([name, '.mean'])];
    all_names = [all_names, ['mean(', eval([name, '.name']),')']];
    all_vars = [all_vars eval([name, '.std'])];
    all_names = [all_names, ['std(', eval([name, '.name']),')']];
    all_vars = [all_vars eval([name, '.skewness'])];
    all_names = [all_names, ['skewness(', eval([name, '.name']),')']];
    all_vars = [all_vars eval([name, '.kurtosis'])];
    all_names = [all_names, ['kurtosis(', eval([name, '.name']),')']];
end
all_vars = [all_vars other_stats];
all_names = [all_names, ...
    'TB(whole)', 'TB(part)', 'Vd(whole)/mean(v)(whole)', ...
    'Vd(part)/mean(v)(whole)', 'Vd(part)/mean(v)(part)', ...
    'tRp/tRm(whole)', 'mean run speed', ...
    'tp/tm', 'dp/dm', 'switching frequency', ...
    'mean(tpm)', 'std(tpm)'];

all_Series = [dfdt, dxdt, tb1, tb21, tb41, tb81];
all_Joint = [dfdt_x_tb1, dfdt_x_tb21, dfdt_x_tb41, dfdt_x_tb81, dxdt_x_tb1, dxdt_x_tb21, dxdt_x_tb41, dxdt_x_tb81];

n_input = 2; % number of input Series
n_output = 4; % number of output Series
s1 = 1; % index of initial input Series
s2 = n_input+1; % index of initial output Series

stats_interactiveTrackDisplay(all_vars, all_names, ...
    (s1-1)*(N_PC+4) + N_PC + [1 2 3], (s2-1)*(N_PC+4) + N_PC + [1 2 3], ...
    N_series*(N_PC+4)+ 3, all_Series, [s1 s2], all_Joint, (n_input-1)*n_output+s2-s1, ...
    pos, dt, [num2str(traj_length),'_'])