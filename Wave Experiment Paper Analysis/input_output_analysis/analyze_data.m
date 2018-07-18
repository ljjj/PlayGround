traj_length = 120;
refilter = false;
%% analyze experimental data
% data located in dropbox
% Dropbox (emonetlab)/users/adam_waite/Documenteed
% all_tracks = cell(1,2);s/Projects/gradient/experiments/microfluidics/170222-RP437_gradient_4X/170222-RP437_gradient_4X_001.swimtracker.mat
% Documents/Projects/gradient/experiments/microfluidics/no_gate_4X_all/all_tracks.swimtracker.mat

% initialize others
other_stats = [];
dt = 0.05;
N_times = round(traj_length/dt); % each cut trajectory length
pos = zeros(0, N_times, 2);
c = 0;

% initialize Series and Joint structs
N_bins = 50;
N_lag = 100;
dfdt   = SeriesInitialization('df/dt',        N_bins, -0.01, 0.01, dt, N_times, N_lag);
dxdt   = SeriesInitialization('dx/dt',        N_bins,   -40,   40, dt, N_times, N_lag);
tb1    = SeriesInitialization('tb1',          N_bins,     0,    1, dt, N_times, N_lag);
tb21   = SeriesInitialization('tb21',         N_bins,     0,    1, dt, N_times, N_lag);
tb41   = SeriesInitialization('tb41',         N_bins,     0,    1, dt, N_times, N_lag);
tb81   = SeriesInitialization('tb81',         N_bins,     0,    1, dt, N_times, N_lag);
dfdt_x_tb1    = JointInitialization(dfdt, tb1,    dt, N_lag, N_lag, N_lag);
dfdt_x_tb21   = JointInitialization(dfdt, tb21,   dt, N_lag, N_lag, N_lag);
dfdt_x_tb41   = JointInitialization(dfdt, tb41,   dt, N_lag, N_lag, N_lag);
dfdt_x_tb81   = JointInitialization(dfdt, tb81,   dt, N_lag, N_lag, N_lag);
dxdt_x_tb1    = JointInitialization(dxdt, tb1,    dt, N_lag, N_lag, N_lag);
dxdt_x_tb21   = JointInitialization(dxdt, tb21,   dt, N_lag, N_lag, N_lag);
dxdt_x_tb41   = JointInitialization(dxdt, tb41,   dt, N_lag, N_lag, N_lag);
dxdt_x_tb81   = JointInitialization(dxdt, tb81,   dt, N_lag, N_lag, N_lag);

% load data
data_file = 'C:\scratch\170222-RP437_gradient_4X_001.swimtracker.mat';
gradient_file = 'C:\scratch\170222-gradient.mat';
tstr = '170222-RP437_gradient_4X';
if ~exist('data','var') || isempty(data)
    disp(['loading data ', data_file, ' ...'])
    load(data_file)
    load(gradient_file)
    data.tracks = tracks;
    clear tracks
end

% filter data
if refilter || ~exist('filtered_data','var') || isempty(filtered_data)
    disp('filtering data ...')
    clear filter
    filter.trajtime = [traj_length, Inf]; % filter for long track length
    filter.median_sigma = [0 1.65];
    [filtered_data,ids] = filterData(data, [], filter); % filterData.m is in the repository "Microscope"
end

ratios = runDirectionRatio(filtered_data);
for i = 1:length(filtered_data.tracks)
    track = filtered_data.tracks(i);
    total_ratio = ratios(i);
    % calculated window-averaged tb
    tb1_totalseries    = track.tumble;
    tb1_totalseries(tb1_totalseries == 2) = 0; % remove intermediate state
    tb1_totalseries(isnan(tb1_totalseries)) = 0; % remove nan
    tb21_totalseries   = smooth(tb1_totalseries, 21);
    tb41_totalseries   = smooth(tb1_totalseries, 41);
    tb81_totalseries   = smooth(tb1_totalseries, 81);
    % cut the trajectory
    nc_j = 0; % index offset for each cut trajectory
    while nc_j+N_times+1 <= numel(track.dx) % find pdf for each cut trajectory
        % find df/dt
        dfdt_timeseries = diff(log(gradient(1).F_smoothed(track.x((nc_j+1):(nc_j+N_times+1)), ...
                                                   track.time((nc_j+1):(nc_j+N_times+1)))))/dt;
        dfdt = SeriesAppend(dfdt, dfdt_timeseries);
        % find dx/dt
        dxdt_timeseries = track.dx((nc_j+1):(nc_j+N_times));
        dxdt = SeriesAppend(dxdt, dxdt_timeseries);
        % find tb
        tb1_timeseries   = tb1_totalseries ((nc_j+1):(nc_j+N_times));
        tb21_timeseries  = tb21_totalseries((nc_j+1):(nc_j+N_times));
        tb41_timeseries  = tb41_totalseries((nc_j+1):(nc_j+N_times));
        tb81_timeseries  = tb81_totalseries((nc_j+1):(nc_j+N_times));
        tb1   = SeriesAppend(tb1,  tb1_timeseries);
        tb21  = SeriesAppend(tb21, tb21_timeseries);
        tb41  = SeriesAppend(tb41, tb41_timeseries);
        tb81  = SeriesAppend(tb81, tb81_timeseries);
        % joint calculations
        dfdt_x_tb1    = JointAppend(dfdt_x_tb1,    dfdt, tb1);
        dfdt_x_tb21   = JointAppend(dfdt_x_tb21,   dfdt, tb21);
        dfdt_x_tb41   = JointAppend(dfdt_x_tb41,   dfdt, tb41);
        dfdt_x_tb81   = JointAppend(dfdt_x_tb81,   dfdt, tb81);
        dxdt_x_tb1    = JointAppend(dxdt_x_tb1,    dxdt, tb1);
        dxdt_x_tb21   = JointAppend(dxdt_x_tb21,   dxdt, tb21);
        dxdt_x_tb41   = JointAppend(dxdt_x_tb41,   dxdt, tb41);
        dxdt_x_tb81   = JointAppend(dxdt_x_tb81,   dxdt, tb81);
        % record other stats
        spd = track.speed((nc_j+1):(nc_j+N_times));
        tb = track.tumble((nc_j+1):(nc_j+N_times));
        mrs = nanmean(spd(tb == 0)); % mean run speed
        tb(tb == 2) = 0;
        other_stats = [other_stats; track.tumblebias, nanmean(tb), ... % tumble bias of the whole trajectory or the cut trajectory
            (track.x_raw(end)-track.x_raw(1))/track.trajtime/track.meanrunspeed, ... % normalized drift speed of the whole trajectory
            (track.x_raw(nc_j+N_times)-track.x_raw(nc_j+1))/traj_length/track.meanrunspeed, ... % normalized drift speed of the cut trajectory by the whole trajectory mean run speed
            (track.x_raw(nc_j+N_times)-track.x_raw(nc_j+1))/traj_length/mrs, ... % normalized drift speed of the cut trajectory
            total_ratio, ... % ratio of run up length and run down length of the whole trajectory
            mrs]; % mean run speed
        % record positions
        c = c + 1;
        pos(c,:,:) = [track.x_raw((nc_j+1):(nc_j+N_times)) track.y_raw((nc_j+1):(nc_j+N_times))];
        nc_j = nc_j + N_times; % advance to the next cut trajectory
    end
end

dfdt   = SeriesSummarize(dfdt);
dxdt   = SeriesSummarize(dxdt);
tb1    = SeriesSummarize(tb1);
tb21   = SeriesSummarize(tb21);
tb41   = SeriesSummarize(tb41);
tb81   = SeriesSummarize(tb81);

%% calculate direcitonal bias quantities
time = dt*(0:(N_times-2));
dpos = diff(pos,1,2); % change in position over time
theta = squeeze(atan2(sqrt(sum(dpos(:,:,2:end).^2,3)), dpos(:,:,1) )); % angle made with gradient
posDir = theta <= pi/2; % logical of whether in + direction
dwellTimes = nan(size(posDir)); % dwell times in each direction
changeState = diff(posDir,1,2); % indicate where a change in state happens
dWTcounts = zeros(2,size(posDir,2)); % wait time distribution, first index is + / - direction
displacements = zeros(2,floor(size(posDir,1)*size(posDir,2)/2)); % displacement values, first index is + / - direction
dInd = [1 1]; % index to add displacement values
DB = mean(posDir, 2); % direction bias (fraction of time spent in positive directions)
DBdisp = zeros(size(posDir,1),1); % direction bias of displacement (fraction of distance travelled in positive directions)
SF = zeros(size(posDir,1),1); % switching frequency between + / - directions
for n = 1:size(posDir,1)
    changes = changeState(n,:);
    changeInd = find(changes ~= 0);
    totalDis = [0 0];
    for k = 2:numel(changeInd)
        state = changes(changeInd(k)) == 1; % is the change from - to +
        len = changeInd(k) - changeInd(k-1); % time of staying in this direction
        dwellTimes(n, (changeInd(k-1)+1):(changeInd(k)+1)) = - len * (state-0.5)*2;
        dWTcounts(state+1, len) = dWTcounts(state+1, len) + 1;
        displace = pos(n,changeInd(k)+1,1) - pos(n,changeInd(k-1)+1,1); % displacement in this direction
        displacements(state+1, dInd(state+1)) = displace;
        totalDis(state+1) = totalDis(state+1) + abs(displace);
        dInd(state+1) = dInd(state+1) + 1;
    end
    DBdisp(n) = totalDis(1) / sum(totalDis);
    SF(n) = numel(changeInd)/time(end);
end
displacements(:, (max(dInd)+1):end) = [];

other_stats = [other_stats DB DBdisp SF ...
    nanmean(dwellTimes, 2), nanstd(dwellTimes, 0, 2)];
%% save analysis results
save_name = ['stats_',num2str(traj_length),'.mat'];
disp(['saving analysis result ', save_name, ' ...'])
save(save_name, 'traj_length', 'dfdt', 'dxdt', 'tb1', 'tb21', 'tb41', 'tb81', ...
    'dfdt_x_tb1', 'dfdt_x_tb21', 'dfdt_x_tb41', 'dfdt_x_tb81', ...
    'dxdt_x_tb1', 'dxdt_x_tb21', 'dxdt_x_tb41', 'dxdt_x_tb81', ...
    'other_stats','pos','dt','-v7.3')