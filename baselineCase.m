numCrops = 1000; 
initialI = 0; initialC = 1; 
r=0.05; b = r/numCrops;
beta = b; 
epsilon = 0.015; 
gamma = 452;

numSims = 1000; tFinal = 5000; progress = "yes"; numRuns = 1000; newFigure = "yes";

popSize = numCrops;

simData = runSpread_SCI_1(popSize,initialI,initialC,beta,epsilon,gamma,numSims,tFinal,progress);

 
sampleSizeVec = 20:20:100;
sampleIntervalVec = 30:30:120;

EDIstore = zeros(length(sampleSizeVec),length(sampleIntervalVec));
perc95store = zeros(length(sampleSizeVec),length(sampleIntervalVec));

for i=1:length(sampleSizeVec)
    sampleSize = sampleSizeVec(i)
    for j=1:length(sampleIntervalVec)
        sampleInterval = sampleIntervalVec(j);
        
        [sampleData, EDI, EDT, EDR, perc95] = runSampling_SCI_1(simData,numRuns,popSize,sampleSize,sampleInterval,tFinal,progress);
        
        EDIstore(i,j) = EDI;
        perc95store(i,j) = perc95;
        
    end
end

%%
figure();
% levels = [0 0.5 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40];
[C,h] = contourf(EDIstore'/10);

h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

colbar = colorbar;
ylabel(colbar, 'Baseline EDP (% of population)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

xlabel('Sample size (N)')
ylabel('Sample interval (\Delta days)')

xticks = 1:1:length(sampleSizeVec);
xticklabels = strsplit(num2str(sampleSizeVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

yticks = 1:1:length(sampleIntervalVec);
yticklabels = strsplit(num2str(sampleIntervalVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

box off; set(gca,'Fontsize',16,'Linewidth',2);

%%
figure();
% levels = [0 0.5 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40];
[C,h] = contourf(perc95store'/10);

h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

colbar = colorbar;
ylabel(colbar, 'Baseline DP_{95} (% of population)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

xlabel('Sample size (N)')
ylabel('Sample interval (\Delta days)')

xticks = 1:1:length(sampleSizeVec);
xticklabels = strsplit(num2str(sampleSizeVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

yticks = 1:1:length(sampleIntervalVec);
yticklabels = strsplit(num2str(sampleIntervalVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

box off; set(gca,'Fontsize',16,'Linewidth',2);







        