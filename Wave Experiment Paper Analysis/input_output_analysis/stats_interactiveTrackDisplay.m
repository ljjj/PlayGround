function stats_interactiveTrackDisplay(all_vars, all_names, ...
    p1_inds, p2_inds, c_ind, all_Series, s_inds, all_Joint, j_ind, ...
    pos, dt, save_pref)
% This function interactively compares different all_vars and 
% show trajectory plots, Series plots and Joint plots for save.
% all_vars  : array of variables (N_tracks x N_vars) to compare
% all_names : cell array of corresponding names of all_vars
% p1_inds   : array 1x3 of indices of variables to show in panel 1
% p2_inds   : array 1x3 of indices of variables to show in panel 2
% c_ind     : index of variables to show in colorbar
% all_Series: array of structs defined in SeriesInitialization.m
% s_inds    : array 1x2 of indices of Series structs to show
% all_Joint : array of structs defined in JointInitialization.m
% j_ind     : index of Joint struct to show
% pos       : N_tracks x N_times x 3(or 2) array of positions of track observations
% dt        : frame interval
% save_pref : name prefix to save the plots

% input consistency check
[N_tracks, N_vars] = size(all_vars);
N_Series = numel(all_Series);
N_Joint = numel(all_Joint);
if numel(all_names) ~= N_vars
    error('variable names mismatch')
end
for s = 1:N_Series
    if any(size(all_Series(s).n) ~= [N_tracks, all_Series(s).N_bins])
        error(['Series ', num2str(s),' size mismatch'])
    end
end
if numel(p1_inds) ~= 3
    error('p1_inds dimension wrong')
end
if numel(p2_inds) ~= 3
    error('p2_inds dimension wrong')
end
if numel(c_ind) ~= 1
    error('c_ind dimension wrong')
end
if numel(s_inds) ~= 2
    error('s_inds dimension wrong')
end
if numel(j_ind) ~= 1
    error('j_ind dimension wrong')
end
if size(pos,1) ~= N_tracks
    error('position size wrong')
end
N_times = size(pos,2);

% initialization
global scatttered1 scatttered2 % scatter plot handlers
global c1Str c2Str % scatter plot choices
global xMax xMin yMax yMin % positions of the highlighted box
global highlightedBox highlightedPoints boxSelect % highlight handlers
global track_inds % indices in N_tracks of the highlighted points
global track_inds_show % indices of tracks shown in Series and Joint plots
global N_show % number of tracks shown in Series and Joint plots
scatttered1 = [];
scatttered2 = [];
c1Str = '.';
c2Str = '.';
c3Str = '.';
xMax = [];
xMin = [];
yMax = [];
yMin = [];
highlightedBox = [];
highlightedPoints = [];
boxSelect = [];
track_inds = [];
track_inds_show = [];
N_show = 0;
pos = pos - repmat(pos(:,1,:), 1, N_times, 1); % rezero positions

figure('Position',[20 15 1500 760])
m = 3; % subplot rows
n = 4; % subplot cols

% UI positions
UITopY = 730;
UITopCol1X = 190;
UITopCol2X = 540;
UITopCol4X = (UITopCol2X-UITopCol1X)*3+UITopCol1X;
UILeftCol1X = 80;
UILeftCol2X = UILeftCol1X+UITopCol2X-UITopCol1X;
UIRow1Y = 510;
UIRow2Y = 400;
UIRow3Y = 150;

% UI control of panel 1
uicontrol('Style', 'popup', ...
    'String', all_names, 'Value', p1_inds(1), ...
    'Position', [UILeftCol1X UIRow1Y 60 20], ...
    'Callback', @reuseDataX1);
uicontrol('Style', 'popup', ...
    'String', all_names, 'Value', p1_inds(2), ...
    'Position', [UILeftCol1X UIRow1Y-25 60 20], ...
    'Callback', @reuseDataY1);
uicontrol('Style', 'popup', ...
    'String', all_names, 'Value', p1_inds(3), ...
    'Position', [UILeftCol1X UIRow1Y-50 60 20], ...
    'Callback', @reuseDataZ1);
xlogCheck1 = uicontrol('Style', 'checkbox', ...
    'String', 'xlog', 'Value', 0,...
    'Position', [UILeftCol1X+65 UIRow1Y 40 20], ...
    'Callback', @checkXLog1);
ylogCheck1 = uicontrol('Style', 'checkbox', ...
    'String', 'ylog', 'Value', 0,...
    'Position', [UILeftCol1X+65 UIRow1Y-25 40 20], ...
    'Callback', @checkYLog1);
zlogCheck1 = uicontrol('Style', 'checkbox', ...
    'String', 'zlog', 'Value', 0,...
    'Position', [UILeftCol1X+65 UIRow1Y-50 40 20], ...
    'Callback', @checkZLog1);

    function reuseDataX1(source, callbackdata)
        p1_inds(1) = source.Value;
        scatter1
    end

    function reuseDataY1(source, callbackdata)
        p1_inds(2) = source.Value;
        scatter1
    end

    function reuseDataZ1(source, callbackdata)
        p1_inds(3) = source.Value;
        scatter1
    end

    function checkXLog1(source, callbackdata)
        subplot(m,n,1)
        if source.Value
            set(gca,'XScale','log')
        else
            set(gca,'XScale','linear')
        end
    end

    function checkYLog1(source, callbackdata)
        subplot(m,n,1)
        if source.Value
            set(gca,'YScale','log')
        else
            set(gca,'YScale','linear')
        end
    end

    function checkZLog1(source, callbackdata)
        subplot(m,n,1)
        if source.Value
            set(gca,'ZScale','log')
        else
            set(gca,'ZScale','linear')
        end
    end

% UI control of panel 2
uicontrol('Style', 'popup', ...
    'String', all_names, 'Value', p2_inds(1), ...
    'Position', [UILeftCol2X UIRow1Y 60 20], ...
    'Callback', @reuseDataX2);
uicontrol('Style', 'popup', ...
    'String', all_names, 'Value', p2_inds(2), ...
    'Position', [UILeftCol2X UIRow1Y-25 60 20], ...
    'Callback', @reuseDataY2);
uicontrol('Style', 'popup', ...
    'String', all_names, 'Value', p2_inds(3), ...
    'Position', [UILeftCol2X UIRow1Y-50 60 20], ...
    'Callback', @reuseDataZ2);
xlogCheck2 = uicontrol('Style', 'checkbox', ...
    'String', 'xlog', 'Value', 0,...
    'Position', [UILeftCol2X+65 UIRow1Y 40 20], ...
    'Callback', @checkXLog2);
ylogCheck2 = uicontrol('Style', 'checkbox', ...
    'String', 'ylog', 'Value', 0,...
    'Position', [UILeftCol2X+65 UIRow1Y-25 40 20], ...
    'Callback', @checkYLog2);
zlogCheck2 = uicontrol('Style', 'checkbox', ...
    'String', 'zlog', 'Value', 0,...
    'Position', [UILeftCol2X+65 UIRow1Y-50 40 20], ...
    'Callback', @checkZLog2);

    function reuseDataX2(source, callbackdata)
        p2_inds(1) = source.Value;
        scatter2
    end

    function reuseDataY2(source, callbackdata)
        p2_inds(2) = source.Value;
        scatter2
    end

    function reuseDataZ2(source, callbackdata)
        p2_inds(3) = source.Value;
        scatter2
    end

    function checkXLog2(source, callbackdata)
        subplot(m,n,2)
        if source.Value
            set(gca,'XScale','log')
        else
            set(gca,'XScale','linear')
        end
    end

    function checkYLog2(source, callbackdata)
        subplot(m,n,2)
        if source.Value
            set(gca,'YScale','log')
        else
            set(gca,'YScale','linear')
        end
    end

    function checkZLog2(source, callbackdata)
        subplot(m,n,2)
        if source.Value
            set(gca,'ZScale','log')
        else
            set(gca,'ZScale','linear')
        end
    end

% UI control of scatter type
scst1 = uicontrol('Style', 'checkbox', ...
    'String', 'circles', ...
    'Position', [UITopCol1X+50 UITopY 60 20], ...
    'Callback', @circle1);
scst2 = uicontrol('Style', 'checkbox', ...
    'String', 'circles', ...
    'Position', [UITopCol2X+50 UITopY 60 20], ...
    'Callback', @circle2);
scst3 = uicontrol('Style', 'checkbox', ...
    'String', 'circles', ...
    'Position', [UITopCol4X+50 UITopY 60 20], ...
    'Callback', @circle3);

    function circle1(source, callbackdata)
        if scst1.Value
            c1Str = 'o';
        else
            c1Str = '.';
        end
        scatter1
    end

    function circle2(source, callbackdata)
        if scst2.Value
            c2Str = 'o';
        else
            c2Str = '.';
        end
        scatter2
    end

    function circle3(source, callbackdata)
        if scst3.Value
            c3Str = 'o';
        else
            c3Str = '.';
        end
        updateSC(all_Series(s_inds(1)), all_Series(s_inds(2)), track_inds_show, 4)
        textAnnotation(N_show)
    end

% UI control of colorbar
uicontrol('Style', 'popup', ...
    'String', all_names, 'Value', c_ind,...
    'Position', [(UITopCol1X+UITopCol2X)/2-20 UITopY 60 20], ...
    'Callback', @relabel);

    function relabel(source, callbackdata)
        c_ind = source.Value;
        scatter1
        scatter2
    end

% UI control of box selection
st = uicontrol('Style', 'popup', ...
    'String', {'manual box', 'auto square'}, ...
    'Position', [(UITopCol1X+UITopCol2X)/2+45 UITopY 60 20]);
uicontrol('Style', 'pushbutton', 'String', 'Select1', ...
    'Position', [UITopCol1X UITopY 45 20], ...
    'Callback', @select1);
uicontrol('Style', 'pushbutton', 'String', 'Select2', ...
    'Position', [UITopCol2X UITopY 45 20], ...
    'Callback', @select2);

    function select1(source, callbackdata)
        subplot(m,n,1)
        selectBox
        track_inds = find(all_vars(:,p1_inds(1)) <= xMax & ...
            all_vars(:,p1_inds(1)) >= xMin & ...
            all_vars(:,p1_inds(2)) <= yMax & ...
            all_vars(:,p1_inds(2)) >= yMin); % obtain track indices
        subplot(m,n,2)
        highlightPoints(all_vars(track_inds,p2_inds(1)), ...
            all_vars(track_inds,p2_inds(2)), ...
            all_vars(track_inds,p2_inds(3)))
        updateTracks
        boxSelect = 1;
        if strcmp(sv.Enable,'off')
            set(sv,'Enable','on');
        end
    end

    function select2(source, callbackdata)
        subplot(m,n,2)
        selectBox
        track_inds = find(all_vars(:,p2_inds(1)) <= xMax & ...
            all_vars(:,p2_inds(1)) >= xMin & ...
            all_vars(:,p2_inds(2)) <= yMax & ...
            all_vars(:,p2_inds(2)) >= yMin); % obtain track indices
        subplot(m,n,1)
        highlightPoints(all_vars(track_inds,p1_inds(1)), ...
            all_vars(track_inds,p1_inds(2)), ...
            all_vars(track_inds,p1_inds(3)))
        updateTracks
        boxSelect = 2;
        if strcmp(sv.Enable,'off')
            set(sv,'Enable','on');
        end
    end

% select and highlight
    function selectBox()
        view(2)
        switch st.Value
            case 1
                % get the box
                [x1, y1] = ginput(1);
                delete(highlightedBox)
                hold on
                cross = plot(x1,y1,'r+','MarkerSize',20); % temporary cross
                [x2, y2] = ginput(1);
                delete(cross)
                xMin = min([x1 x2]);
                xMax = max([x1 x2]);
                yMin = min([y1 y2]);
                yMax = max([y1 y2]);
            case 2
                [x, y] = ginput(1);
                delete(highlightedBox)
                hold on
                boxSize = sqrt(2);
                xMin = x/boxSize;
                xMax = x*boxSize;
                yMin = y/boxSize;
                yMax = y*boxSize;
        end
        highlightedBox = plot([xMin xMax xMax xMin xMin], [yMin yMin yMax yMax yMin], 'r','LineWidth',2);
        hold off
    end

    function highlightPoints(x,y,z)
        delete(highlightedPoints)
        hold on
        highlightedPoints = plot3(x,y,z,'xk','MarkerSize',6);
        hold off
    end

% UI control of tracks to show
tracks_ns = [5 10 20 50 100 200 500 1000];
tracks_strs = {};
for i = 1:numel(tracks_ns)
    tracks_strs = [tracks_strs num2str(tracks_ns(i))];
end
ts = uicontrol('Style', 'popup', ...
    'String', tracks_strs, 'Value', 2, ...
    'Position', [UILeftCol1X UIRow3Y 50 20], ...
    'Callback', @updateTracks);

    function updateTracks(source, callbackdata)
        N_show = tracks_ns(ts.Value); % get the number of tracks to show
        if numel(track_inds) > N_show
            track_inds_show = track_inds(1:N_show); % remove tracks more than enough
        else
            track_inds_show = track_inds;
        end
        N_show = numel(track_inds_show); % count actual tracks to show
        
        % plot the trajectories in x and y
%         subplot(m,n,3)
%         delete(gca)
%         subplot(m,n,3)
%         if ~isempty(track_inds_show)
%             if size(pos,3) == 3
%                 plot(pos(track_inds_show,:,1)', sqrt(sum(pos(track_inds_show,:,2:end).^2, 3))')
%             else
%                 plot(pos(track_inds_show,:,1)', pos(track_inds_show,:,2)')
%             end
%         end
%         xlabel('x (\mum)')
%         ylabel('y (\mum)')
%         axis equal
%         title('track in space')
%         textAnnotation(N_show)
        
        % plot the trajectories in x vs. t
        subplot(m,n,3)
        delete(gca)
        subplot(m,n,3)
        if ~isempty(track_inds_show)
            x = squeeze(pos(track_inds_show,:,1))';
            t = (0:(N_times-1))*dt;
            plotWithMean(t, x)
        end
        xlabel('time (s)')
        ylabel('x (\mum)')
        axis tight
        title('trajectory in time')
        textAnnotation(N_show)
        
        updateSeries
        updateJoint
    end

    function plotWithMean(x, y, LineSpec)
        if numel(y) == 0
            return
        end
        if ~exist('LineSpec','var')
            LineSpec = '-';
        end
        plot(x, y, LineSpec)
        if size(y, 2) > 1
            hold on
            meanY = mean(y,2);
            stdY = std(y,0,2);
            plot(x, meanY, 'k', 'LineWidth', 2)
            plot(x, meanY + stdY, 'k--', 'LineWidth', 1.4)
            plot(x, meanY - stdY, 'k--', 'LineWidth', 1.4)
            hold off
        end
    end

    function textAnnotation(N_t)
        text(0.1, 0.9, [num2str(N_t), ' tracks'], 'Units', 'normalized')
    end

% UI control of Series and Joint to show
Series_names = {};
for s = 1:N_Series
    Series_names = [Series_names, all_Series(s).name];
end
Joint_names = {};
for s = 1:N_Joint
    Joint_names = [Joint_names, [all_Joint(s).name1, ' x ', all_Joint(s).name2]];
end
s1 = uicontrol('Style', 'popup', ...
    'String', Series_names, 'Value', s_inds(1), ...
    'Position', [UILeftCol1X UIRow3Y+75 50 20], ...
    'Callback', @updateSeries);
s2 = uicontrol('Style', 'popup', ...
    'String', Series_names, 'Value', s_inds(2), ...
    'Position', [UILeftCol1X UIRow3Y+50 50 20], ...
    'Callback', @updateSeries);
jt = uicontrol('Style', 'popup', ...
    'String', Joint_names, 'Value', j_ind, ...
    'Position', [UILeftCol1X UIRow3Y+25 50 20], ...
    'Callback', @updateJoint);

    function updateSeries(source, callbackdata)
        s_inds(1) = s1.Value;
        s_inds(2) = s2.Value;
        
        % plot time series
%         updateTimeSeries(all_Series(s_inds(1)),track_inds_show,5)
%         textAnnotation(N_show)
%         updateTimeSeries(all_Series(s_inds(2)),track_inds_show,6)
%         textAnnotation(N_show)
        
        % plot power spectrum
        updatePSD(all_Series(s_inds(1)),track_inds_show,5)
        textAnnotation(N_show)
        updatePSD(all_Series(s_inds(2)),track_inds_show,6)
        textAnnotation(N_show)
        
        % plot pdf
        updatePDF(all_Series(s_inds(1)),track_inds_show,7)
        textAnnotation(N_show)
        updatePDF(all_Series(s_inds(2)),track_inds_show,8)
        textAnnotation(N_show)
        
        % plot acf and pacf
        updateACFPACF(all_Series(s_inds(1)),track_inds_show,9)
        textAnnotation(N_show)
        updateACFPACF(all_Series(s_inds(2)),track_inds_show,10)
        textAnnotation(N_show)
        
        % plot scatter
        updateSC(all_Series(s_inds(1)), all_Series(s_inds(2)), track_inds_show, 4)
        textAnnotation(N_show)
    end

    function updateJoint(source, callbackdata)
        j_ind = jt.Value;
        
        % plot cross-correlation
        updateXC(all_Joint(j_ind), track_inds_show, 11)
        textAnnotation(N_show)
        
        % plot mutual information
        updateMI(all_Joint(j_ind), track_inds_show, 12)
        textAnnotation(N_show)
        
        % plot linear filter
%         updateLF(all_Joint(j_ind), track_inds_show, 4)
%         textAnnotation(N_show)
    end

% update Series plots
    function updateTimeSeries(Series, inds, k)
        subplot(m,n,k)
        delete(gca)
        subplot(m,n,k)
%         plotWithMean(Series.time, Series.timeseries(inds,:)')
        plot(Series.time, Series.timeseries(inds,:)')
        xlabel('time (s)')
        ylabel(Series.name)
        axis tight
    end

    function updatePSD(Series, inds, k)
        subplot(m,n,k)
        delete(gca)
        subplot(m,n,k)
        hold on
%         plotWithMean(Series.freq, Series.psd(inds,:)')
        for i = inds(:)'
            plot(Series.freq, smooth(Series.psd(i,:), smooth_windows(psd_smooth.Value)))
        end
        plot(Series.freq, 1e-3./Series.freq, '--k')
        plot(Series.freq, 1e-3./Series.freq.^2, '--r')
        text(0.1, 0.1, '1/f', 'Units', 'normalized')
        text(0.1, 0.2, '1/f^2', 'Units', 'normalized', 'color', 'r')
        hold off
        xlabel('frequency (Hz)')
        ylabel(['PSD(',Series.name,')'])
        if xlogPSDCheck.Value
            set(gca, 'XScale', 'log')
        else
            set(gca, 'XScale', 'linear')
        end
        set(gca, 'YScale', 'log')
        axis tight
    end

    function updatePDF(Series, inds, k)
        subplot(m,n,k)
        delete(gca)
        subplot(m,n,k)
%         plotWithMean(Series.ctrs, Series.n(inds,:)')
        if ~isempty(inds)
            plot(Series.ctrs, Series.n(inds,:)')
        end
        xlabel(Series.name)
        ylabel(['PDF(',Series.name,')'])
        axis tight
    end

    function updateACFPACF(Series, inds, k)
        subplot(m,n,k)
        delete(gca)
        subplot(m,n,k)
        hold on
%         plotWithMean(Series.acf_lags, Series.acf(inds,:)')
%         plotWithMean(Series.acf_lags, Series.pacf(inds,:)', '--')
        if ~isempty(inds)
            plot(Series.acf_lags, Series.acf(inds,:)')
            plot(Series.acf_lags, Series.pacf(inds,:)', '--')
        end
        hold off
        xlabel('lag (s)')
        ylabel(['- ACF(',Series.name,') -- PACF(',Series.name,')'])
        axis tight
        ylim([-1 1])
    end

    function updateSC(Series1, Series2, inds, k)
        subplot(m,n,k)
        delete(gca)
        subplot(m,n,k)
        hold on
        if ~isempty(inds)
            for j = 1:numel(inds)
                scatter(Series1.timeseries(inds(j),:), Series2.timeseries(inds(j),:), c3Str)
            end
        end
        xlabel(Series1.name)
        ylabel(Series2.name)
        axis tight
        if startsWith(Series1.name, 'tb')
            xlim([0 1])
        end
        if startsWith(Series2.name, 'tb')
            ylim([0 1])
        end
    end

% update Joint plots
    function updateXC(Joint, inds, k)
        subplot(m,n,k)
        delete(gca)
        subplot(m,n,k)
%         plotWithMean(Joint.xc_lags, Joint.xc(inds,:)')
        if ~isempty(inds)
            plot(Joint.xc_lags, Joint.xc(inds,:)')
        end
        xlabel('lag (s)')
        ylabel(['xcorr(',Joint.name1,',',Joint.name2,')'])
        axis tight
        ylim([-1 1])
    end

    function updateMI(Joint, inds, k)
        subplot(m,n,k)
        delete(gca)
        subplot(m,n,k)
%         plotWithMean(Joint.mi_lags, Joint.mi(inds,:)')
        if ~isempty(inds)
            plot(Joint.mi_lags, Joint.mi(inds,:)')
        end
        xlabel('lag (s)')
        ylabel(['MI(',Joint.name1,',',Joint.name2,')'])
        axis tight
        yl = get(gca, 'YLim');
        yl(1) = 0;
        ylim(yl)
    end

    function updateLF(Joint, inds, k)
        subplot(m,n,k)
        delete(gca)
        subplot(m,n,k)
%         plotWithMean(Joint.lf_lags, Joint.lf(inds,:)')
        if ~isempty(inds)
            plot(Joint.lf_lags, Joint.lf(inds,:)')
%             colorOrder = get(gca, 'ColorOrder');
%             N_colors = size(colorOrder,1);
%             for j = i:numel(inds)
%                 text(num2str(Joint.lf_r2(j)), colorOrder(mod(j, N_colors),:))
%             end
        end
        xlabel('lag (s)')
        ylabel('LinearFilter')
        axis tight
        
        axes('Position',[.88 .88 .07 .085])
        box on
        if ~isempty(inds)
            plot(repmat(Joint.lf_r2(inds)',2,1), repmat([0;1],1,numel(inds)))
        end
        xlabel('r^2')
        axis tight
    end

% UI control to smooth PSD
smooth_windows = [1 3 5 9 15 25];
smooth_strs = {};
for i = 1:numel(smooth_windows)
    smooth_strs = [smooth_strs num2str(smooth_windows(i))];
end
psd_smooth = uicontrol('Style', 'popup', ...
    'String', smooth_strs, 'Value', 1, ...
    'Position', [UILeftCol1X UIRow2Y 50 20], ...
    'Callback', @updateSeries);
xlogPSDCheck = uicontrol('Style', 'checkbox', ...
    'String', 'xlog', 'Value', 1,...
    'Position', [UILeftCol1X UIRow2Y-25 40 20], ...
    'Callback', @checkXLogPSD);

    function checkXLogPSD(source, callbackdata)
        subplot(m,n,5)
        if source.Value
            set(gca,'XScale','log')
        else
            set(gca,'XScale','linear')
        end
        subplot(m,n,6)
        if source.Value
            set(gca,'XScale','log')
        else
            set(gca,'XScale','linear')
        end
    end

% UI control of saving plot
sv = uicontrol('Style', 'pushbutton', 'String', 'Save', ...
    'Position', [UILeftCol1X UITopY 40 20], ...
    'Callback', @saveplot, 'Enable','off');

    function saveplot(source, callbackdata)
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 80 90]);
        selectStr = ['[',num2str(xMin),'_',num2str(xMax),'_',num2str(yMin),'_',num2str(yMax),']'];
        saveName = [all_Series(s_inds(1)).name,'_',all_Series(s_inds(2)).name,'_', ...
            all_Joint{j_ind}.name1, 'x', all_Joint{j_ind}.name2, '_', ...
            save_pref,'interactiveTracks_', ...
            all_names{p1_inds(1)},'_',all_names{p1_inds(2)},'_',all_names{p1_inds(3)},'_', ...
            all_names{p2_inds(1)},'_',all_names{p2_inds(2)},'_',all_names{p2_inds(3)},'_', ...
            all_names{c_ind},'_',num2str(boxSelect),'_',selectStr,'_1',c1Str,'_2',c2Str,'.png'];
        saveName = strrep(saveName,'/','');
        disp(['Saving file ',saveName, ' ...'])
        saveas(gcf, saveName)
    end

% show plots
scatter1
scatter2

    function scatter1()
        subplot(m,n,1)
        delete(scatttered1)
        delete(highlightedBox)
        delete(highlightedPoints)
        hold on
        scatttered1 = scatter3(all_vars(:,p1_inds(1)), all_vars(:,p1_inds(2)), all_vars(:,p1_inds(3)), ...
            [], all_vars(:,c_ind), c1Str);
        highlightPoints(all_vars(track_inds,p1_inds(1)), ...
            all_vars(track_inds,p1_inds(2)), ...
            all_vars(track_inds,p1_inds(3)))
        hold off
        xlabel(all_names{p1_inds(1)})
        ylabel(all_names{p1_inds(2)})
        zlabel(all_names{p1_inds(3)})
        cb = colorbar;
        title(cb,all_names{c_ind});
        if xlogCheck1.Value
            set(gca,'XScale','log')
        else
            set(gca,'XScale','linear')
        end
        if ylogCheck1.Value
            set(gca,'YScale','log')
        else
            set(gca,'YScale','linear')
        end
        if zlogCheck1.Value
            set(gca,'ZScale','log')
        else
            set(gca,'ZScale','linear')
        end
        axis tight
    end

    function scatter2()
        subplot(m,n,2)
        delete(scatttered2)
        delete(highlightedBox)
        delete(highlightedPoints)
        hold on
        scatttered2 = scatter3(all_vars(:,p2_inds(1)), all_vars(:,p2_inds(2)), all_vars(:,p2_inds(3)), ...
            [], all_vars(:,c_ind), c2Str);
        highlightPoints(all_vars(track_inds,p2_inds(1)), ...
            all_vars(track_inds,p2_inds(2)), ...
            all_vars(track_inds,p2_inds(3)))
        hold off
        xlabel(all_names{p2_inds(1)})
        ylabel(all_names{p2_inds(2)})
        zlabel(all_names{p2_inds(3)})
        cb = colorbar;
        title(cb,all_names{c_ind});
        if xlogCheck2.Value
            set(gca,'XScale','log')
        else
            set(gca,'XScale','linear')
        end
        if ylogCheck2.Value
            set(gca,'YScale','log')
        else
            set(gca,'YScale','linear')
        end
        if zlogCheck2.Value
            set(gca,'ZScale','log')
        else
            set(gca,'ZScale','linear')
        end
        axis tight
    end
end