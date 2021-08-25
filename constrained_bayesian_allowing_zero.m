% addpath('/Users/francescalovell-read/Documents/DPhil/Sentinels project/Canary_Lab','-end');
%%

numCrops = 1000; 
initialI = 0; initialC = 1; 
r=0.05; b = r/numCrops;
beta = b; betaC = b; betaS = b;
epsilon = 0.015; epsilonC = 0.015; epsilonS = 0.5;
% epsilon = 0.1; epsilonC = 0.1; epsilonS = 0.1;
% epsilon = 1; epsilonC = 1; epsilonS = 1;

gamma = 452; gammaC = 452; gammaS = 49;

sampleSize = 60; sampleInterval = 30;

numSims = 1000; tFinal = 5000; progress = "yes"; numRuns = 1000; newFigure = "yes";

%%
% COMPUTE BASELINE EDI (NO SENTINELS)
simData = runSpread_SCI_1(numCrops,initialI,initialC,beta,epsilon,gamma,20000,tFinal,progress);
[sampleData, EDI, EDT, EDR, perc95] = runSampling_SCI_1(simData,20000,numCrops,sampleSize,sampleInterval,tFinal,progress);
baseline_EDI = EDI;
baseline_95 = perc95;

%%
close all
progress = "no";

% totalSentinels = optimizableVariable('Stot',[0,250],'Type','integer');
totalSentinels = optimizableVariable('Stot',[0,25],'Type','integer');
sampleSent = optimizableVariable('Ssamp',[0,sampleSize],'Type','integer');

% OPTIMISE!

% fun = @(z)testfunc(numCrops+z.Stot,z.Stot,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress,numRuns,sampleSize-z.Ssamp,z.Ssamp,sampleInterval,baseline_EDI);
fun = @(z)testfunc(numCrops+10*z.Stot,10*z.Stot,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress,numRuns,sampleSize-z.Ssamp,z.Ssamp,sampleInterval,baseline_EDI);

results = bayesopt(fun,[totalSentinels,sampleSent],'IsObjectiveDeterministic',false,'MaxObjectiveEvaluations',100,'AcquisitionFunctionName','expected-improvement-plus','XConstraintFcn',@xconstraint,'UseParallel',false);

%%
optimalParams = results.XAtMinEstimatedObjective;
optimalSent = 10*optimalParams.Stot
optimalSamp = optimalParams.Ssamp
valueAtOpt = results.MinEstimatedObjective

new_EDI = 0.01*valueAtOpt*baseline_EDI+baseline_EDI;

params = [sampleSize, sampleInterval, optimalSent, optimalSamp, valueAtOpt, baseline_EDI, new_EDI];

% sS = num2str(sampleSize);
% sI = num2str(sampleInterval);
% filename = strcat('sS=',sS,',',32,'sI=',sI);

% filename = 'test3';
% dlmwrite(filename,params,'-append');
% % dlmwrite(filename,optimalSent,'-append');
% % dlmwrite(filename,optimalSamp,'-append');
% % dlmwrite(filename,valueAtOpt,'-append');


%%
% Use test2b for Fig 6
% Use test3b for Fig 7
mymat = dlmread('test2b');
col1 = mymat(:,1);
col2 = mymat(:,2);
col3 = mymat(:,3);
col4 = mymat(:,4);
col5 = mymat(:,5);
col6 = mymat(:,6);
col7 = mymat(:,7);
col8 = mymat(:,8);

%% PLOT OPTIMAL SENTINEL NUMBER!
optSentNum = reshape(col3,[4 5]);
figure(); 

[C,h] = contourf(optSentNum);
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% caxis([0 100])
colbar = colorbar;
ylabel(colbar, 'Optimal number of sentinels', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

set(gca,'fontsize',16,'Linewidth',2); box off;

xlabel('Sample size (N)')
ylabel('Sample interval (\Delta days)')

xticks = 1:1:5;
xticklabels = strsplit(num2str([20 40 60 80 100]));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

yticks = 1:1:4;
yticklabels = strsplit(num2str([30 60 90 120]));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

%% PLOT OPTIMAL SENTINEL PROPORTION!
optSentProp = reshape(col4 ,[4 5])./[20 40 60 80 100];
figure(); 

% levels = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
[C,h] = contourf(optSentProp);
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

caxis([0 1])

colbar = colorbar;
ylabel(colbar, 'Optimal proportion of sentinels in sample', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

set(gca,'fontsize',16,'Linewidth',2); box off;

xlabel('Sample size (N)')
ylabel('Sample interval (\Delta days)')

xticks = 1:1:5;
xticklabels = strsplit(num2str([20 40 60 80 100]));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

yticks = 1:1:4;
yticklabels = strsplit(num2str([30 60 90 120]));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

%% PLOT EDI REDUCTION!
EDIreduction = reshape(col5 ,[4 5]);
figure(); 

[C,h] = contourf(EDIreduction);
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end
% caxis([0 100])

colbar = colorbar;
ylabel(colbar, 'Change in EDP from baseline at optimum (%)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

set(gca,'fontsize',16,'Linewidth',2); box off;

xlabel('Sample size (N)')
ylabel('Sample interval (\Delta days)')

xticks = 1:1:5;
xticklabels = strsplit(num2str([20 40 60 80 100]));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

yticks = 1:1:4;
yticklabels = strsplit(num2str([30 60 90 120]));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

%% PLOT RESULTANT EDP
EDIreduction = reshape(col8 ,[4 5]);
figure(); 

[C,h] = contourf(EDIreduction);
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end
% caxis([0 100])

colbar = colorbar;
ylabel(colbar, 'EDP at optimum (% of total population)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

set(gca,'fontsize',16,'Linewidth',2); box off;

xlabel('Sample size (N)')
ylabel('Sample interval (\Delta days)')

xticks = 1:1:5;
xticklabels = strsplit(num2str([20 40 60 80 100]));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

yticks = 1:1:4;
yticklabels = strsplit(num2str([30 60 90 120]));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);



%%
% % OPTMATSAVE
% % levels = [-50 -40 -30 -20 -10 -5 -1 0 5 10 20 30 40 50];
% figure();
% [C,h] = contourf(valueAtOptMat');
% % [C,h] = contourf(Mstore'-M2store');
% 
% 
% h.LineColor = [1 0 0]; h.LineStyle = ':'; h.LineWidth = 2;
% mylabels = clabel(C,h,'manual','color','k','FontSize',16);
% for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end
% 
% colbar = colorbar;
% ylabel(colbar, 'Percentage change in EDI from baseline', 'Fontsize', 16);
% set(colbar,'linewidth',2,'fontsize',16);
% 
% xlabel('Number of sentinels added');
% ylabel('Sentinel asymptomatic period');
% 
% xticks = 1:1:length(numSentinelsVec);
% xticklabels = strsplit(num2str(numSentinelsVec));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);
% 
% yticks = 1:1:length(gammaSVec);
% yticklabels = strsplit(num2str(gammaSVec));
% set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);
% 
% box off; set(gca,'Fontsize',16,'Linewidth',2);

%%
% % OPTMATSAVE
% % levels = [0 10 20 30 40 50 60 70 80 90];
% figure();
% [C,h] = contourf((100*optimalSentMat)');
% % [C,h] = contourf(Mstore'-M2store');
% 
% 
% h.LineColor = [1 0 0]; h.LineStyle = ':'; h.LineWidth = 2;
% mylabels = clabel(C,h,'manual','color','k','FontSize',16);
% for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end
% 
% colbar = colorbar;
% ylabel(colbar, 'Percentage change in EDI from baseline', 'Fontsize', 16);
% set(colbar,'linewidth',2,'fontsize',16);
% 
% xlabel('Number of sentinels added');
% ylabel('Sentinel asymptomatic period');
% 
% xticks = 1:1:length(numSentinelsVec);
% xticklabels = strsplit(num2str(numSentinelsVec));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);
% 
% yticks = 1:1:length(gammaSVec);
% yticklabels = strsplit(num2str(gammaSVec));
% set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);
% 
% box off; set(gca,'Fontsize',16,'Linewidth',2);

%% CUBIC INTERPOLATION
% [X,Y] = meshgrid(gammaSVec,numSentinelsVec);
% gVec = 30:5:390;
% sVec = 10:1:250;
% [X2,Y2] = meshgrid(gVec,sVec);
% Vq = interp2(X,Y,valueAtOptMat,X2,Y2,'cubic');
% % Vq = interp2(X,Y,optimalSentMat,X2,Y2,'linear');
% %
% figure()
% [C,h] = contourf(Vq');
% 
% h.LineColor = [1 0 0]; h.LineStyle = ':'; h.LineWidth = 2;
% mylabels = clabel(C,h,'manual','color','k','FontSize',16);
% for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end
% 
% colbar = colorbar;
% ylabel(colbar, 'Percentage change in EDI from baseline', 'Fontsize', 16);
% set(colbar,'linewidth',2,'fontsize',16);
% 
% xlabel('Number of sentinels added');
% ylabel('Sentinel asymptomatic period');
% 
% xticks = linspace(1,length(sVec),length(numSentinelsVec));
% xticklabels = strsplit(num2str(numSentinelsVec));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);
% 
% yticks = linspace(1,length(gVec),length(gammaSVec));
% yticklabels = strsplit(num2str(gammaSVec));
% set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);
% 
% box off; set(gca,'Fontsize',16,'Linewidth',2);
% 
% matMin = min(min(Vq));
% [x,y] = find(Vq==matMin);
% optSent = sVec(x)
% optG = gVec(y)

%%
% Cstore=contourc(Vq',[0 0]);
% 
% % Plot!
% figure(); hold on; box on; 
% % set(gca,'fontsize',16,'ColorOrder',parula(length(gammaVec)+1)); 
% grid on;
% for k=1:1
%     myContour = Cstore;
%     myPlot = plot(myContour(1,2:end),myContour(2,2:end));
%     myPlot.LineWidth = 2;
% %     legEntries{k} = strcat('\gamma=',num2str(gammaVec(k)));
% end
% 
% xlabel('Number of sentinels added');
% ylabel('Sentinel asymptomatic period');
% 
% xticks = linspace(1,length(sVec),length(numSentinelsVec));
% xticklabels = strsplit(num2str(numSentinelsVec));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);
% 
% yticks = linspace(1,length(gVec),length(gammaSVec));
% yticklabels = strsplit(num2str(gammaSVec));
% set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);
% 
% 
% xlim([1,length(sVec)]);
% ylim([1,length(gVec)]);
% 
% 
% % leg = legend(legEntries);
% % leg.Title.String = 'Asymptomatic period (days)';
% % leg.FontSize = 16; leg.Location = 'northwest'; leg.Box = 'off'; 
% % leg.NumColumns = 2; leg.Orientation = 'horizontal'; 


%%
function OUTPUT = testfunc(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress,numRuns,cropSampleSize,sentinelSampleSize,sampleInterval,baseline_EDI)
    if numSentinels==0
        OUTPUT = 0;
    else
    simData = runSpread_SCI_2(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress);

    [~, EDI, ~, ~, ~, ~] = runSampling_SCI_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress);
%     [~, ~, ~, ~, ECDI, ~] = runSampling_SCI_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress);

     OUTPUT = 100*(EDI-baseline_EDI)/baseline_EDI;
%      OUTPUT = 100*(ECDI-baseline_EDI)/baseline_EDI;

    end
end

function tf = xconstraint(X)
tf = X.Ssamp<=10*X.Stot;
% tf = X.Ssamp<=X.Stot;
end
