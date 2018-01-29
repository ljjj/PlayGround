%% load
dropbox_dir = 'C:\Users\jl2345\Dropbox (emonetlab)\';
% dropbox_dir = 'D:\Dropbox (emonetlab)\';
raw_data_file = [dropbox_dir, 'paper_wave_diversity\manuscript\Figures\final figure May2016\fig2\fig2_material_May2017\fig2_rawdata.mat'];
load(raw_data_file)
%% plot
Ne = 9; % number of experiments
e_cls = {'r', [1 0.843137254901961 0], 'b'};
q_cls = {'k', 'm', 'c'};
hold on
plot([0 0.7], 0.3*[1 1], 'LineStyle', '--', 'Color', q_cls{1})
plot([0 0.7], 0.1*[1 1], 'LineStyle', '--', 'Color', q_cls{2})
plot([0 0.7], 0.05*[1 1], 'LineStyle', '--', 'Color', q_cls{3})
text(0.1, 0.325, '30%','Color', q_cls{1},'FontSize',12)
text(0.1, 0.125, '10%','Color', q_cls{2},'FontSize',12)
text(0.1, 0.075, '5%','Color', q_cls{3},'FontSize',12)
h_cdf = cell(1:Ne);
h_pdf = cell(1:Ne);
for e = 1:Ne
    if e == 2
%         continue
    end
    subplot(5,1,2)
    hold on
    h_cdf{e} = histogram(data2(e).TB, 60, 'DisplayStyle', 'stairs',...
        'EdgeColor', e_cls{floor((e-1)/3)+1}, 'BinLimits', [0,1], 'Normalization', 'cdf');
    xlabel('TB','FontSize',12)
    ylabel('CDF','FontSize',12)
    ylim([0 1])
    
    subplot(5,1,3)
    hold on
    h_pdf{e} = histogram(data2(e).TB, 60, 'DisplayStyle', 'stairs',...
        'EdgeColor', e_cls{floor((e-1)/3)+1}, 'BinLimits', [0,1], 'Normalization', 'pdf');
    xlabel('TB','FontSize',12)
    ylabel('PDF','FontSize',12)
end
bins = h_cdf{1}.BinEdges(1:end-1) + h_cdf{1}.BinWidth/2;
cdf_mean = cell(1:3);
cdf_std = cell(1:3);
subplot(5,1,1)
hold on
for e = 1:3
    values = [h_cdf{3*e-2}.Values; h_cdf{3*e-1}.Values; h_cdf{3*e}.Values];
    cdf_mean{e} = 1 - mean(values);
    cdf_std{e} = std(values);
    x = bins;
    y = cdf_mean{e};
    er = cdf_std{e};
    plot(x, y, 'Color', e_cls{e}, 'LineWidth', 2)
    patch([x flip(x)],[y+er flip(y-er)], e_cls{e},'EdgeColor','none','LineWidth',1,'EdgeAlpha',0.3,'FaceAlpha',0.3)
end
set(gca,'FontSize',12);
xlim([0 1])
ylim([0 1])
xlabel('TB','FontSize',12)
ylabel('1-CDF','FontSize',12)
%% calculate "escape rate"
phi_mean = cell(1:3);
phi_std = cell(1:3);
for e = 1:3
    values = [h_pdf{3*e-2}.Values./(1-h_cdf{3*e-2}.Values);...
        h_pdf{3*e-1}.Values./(1-h_cdf{3*e-1}.Values);...
        h_pdf{3*e}.Values./(1-h_cdf{3*e}.Values)];
    phi_mean{e} = mean(values);
    phi_std{e} = std(values);
    x = bins;
    y = phi_mean{e};
    er = phi_std{e};
    
    subplot(5,1,5)
    hold on
    plot(x, y, 'Color', e_cls{e}, 'LineWidth', 2)
    patch([x flip(x)],[y+er flip(y-er)], e_cls{e},'EdgeColor','none','LineWidth',1,'EdgeAlpha',0.3,'FaceAlpha',0.3)
    ylim([0 30])
    xlabel('TB','FontSize',12)
    ylabel('mean PDF/(1-CDF)','FontSize',12)

    subplot(5,1,4)
    hold on
    plot(x, values, 'Color', e_cls{e}, 'LineWidth', 2)
    ylim([0 30])
    xlabel('TB','FontSize',12)
    ylabel('individual PDF/(1-CDF)','FontSize',12)
end