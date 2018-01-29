load('fig2_rawdata.mat')
load('jun_analysis\final\wave_speeds_and_cell_counts.mat')

tumbiasgrad=0:0.025:0.7;
Tb50=[];Tb100=[];Tb200=[];ori=[];
for i = 1:length(data2)
    if ~isempty(data2(i).TBdist)
        tb = data2(i).TB;
        [x,y]=hist(tb(tb<0.7),tumbiasgrad);
        dist = x./sum(x)/mean(diff(y));
        distx= y;
        if data2(i).AspConc ==50
            Tb50=[Tb50;dist];
        end
        if data2(i).AspConc ==100
            Tb100=[Tb100;dist];
        end
        if data2(i).AspConc ==200
            Tb200=[Tb200;dist];
        end
        if data2(i).AspConc ==0
            ori=[ori;dist];
        end
    end
end

Ntot = 12.3*10^4;
condP50 = Tb50./ori.*repmat(total_cells(1:3)'/Ntot, 1, numel(tumbiasgrad));
condP100 = Tb100./ori.*repmat(total_cells(4:6)'/Ntot, 1, numel(tumbiasgrad));
condP200 = Tb200./ori.*repmat(total_cells(7:9)'/Ntot, 1, numel(tumbiasgrad));
figure
hold on
plot(tumbiasgrad, condP50, 'LineStyle', '--', 'Color', 'r')
plot(tumbiasgrad, condP100, 'LineStyle', '--', 'Color', [1 0.843137254901961 0])
plot(tumbiasgrad, condP200, 'LineStyle', '--', 'Color', 'b')
plot(tumbiasgrad, mean(condP50), 'LineWidth', 3, 'Color', 'r')
plot(tumbiasgrad, mean(condP100), 'LineWidth', 3, 'Color', [1 0.843137254901961 0])
plot(tumbiasgrad, mean(condP200), 'LineWidth', 3, 'Color', 'b')
ylim([0 2])
xlabel('Tumble bias')
ylabel('P(In the wave|TB)')
