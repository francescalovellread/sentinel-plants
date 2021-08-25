% addpath('/Users/francescalovell-read/Documents/DPhil/Sentinels project/Canary_Lab','-end');
%%

numCrops = 1000; 
initialI = 0; initialC = 1; 
r=0.05; b = r/numCrops;
beta = b; betaC = b; betaS = b;
epsilon = 0.015; epsilonC = 0.015; epsilonS = 0.1;
% epsilon = 1; epsilonC = 1; epsilonS = 1;

gamma = 452; gammaC = 452; 
sampleSize = 20; sampleInterval = 30;
% sampleSize = 100; sampleInterval = 30;

numSims = 1000; tFinal = 5000; progress = "yes"; numRuns = 1000; newFigure = "yes";


% COMPUTE BASELINE EDI (NO SENTINELS)
simData = runSpread_SCI_1(numCrops,initialI,initialC,beta,epsilon,gamma,20000,tFinal,progress);
[sampleData, EDI, EDT, EDR, perc95] = runSampling_SCI_1(simData,20000,numCrops,sampleSize,sampleInterval,tFinal,progress);
baseline_EDI = EDI
% baseline_95 = perc95

% %%
% incidenceDist(sampleData,EDI,0,numCrops,20,"yes");
% cumulativeDist(sampleData,numCrops,"bars");

%%

numSentinelsVec = 10:80:250; popSizeVec = numCrops+numSentinelsVec;
gammaSVec = 30:120:390;

valueAtOptMat = zeros(length(numSentinelsVec),length(gammaSVec));
optimalSentMat = zeros(length(numSentinelsVec),length(gammaSVec));
optimalSentMat2 = zeros(length(numSentinelsVec),length(gammaSVec));

for j = 1:length(numSentinelsVec)
    numSentinels = numSentinelsVec(j)
    popSize = popSizeVec(j);
    
    for k = 1:length(gammaSVec)
        gammaS = gammaSVec(k)


maxNoSentinels = min(numSentinels,sampleSize);

if maxNoSentinels == 0
%     fprintf('No sentinels. Moving on...\n');
    optimalSentMat(j,k) = 0;
    optimalSentMat2(j,k) = 0;
    valueAtOptMat(j,k) = 0;
else
    
numSent = optimizableVariable('S',[0,maxNoSentinels],'Type','integer');

% OPTIMISE!
simData = runSpread_SCI_2(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,1000,tFinal,progress);

fun = @(z)testfunc(simData,1000,popSize,numSentinels,sampleSize-z.S,z.S,sampleInterval,tFinal,progress,baseline_EDI);

results = bayesopt(fun,numSent,'IsObjectiveDeterministic',false,'MaxObjectiveEvaluations',30,'AcquisitionFunctionName','expected-improvement-plus','PlotFcn',[]);

optimalParams = results.XAtMinEstimatedObjective;
optimalSent = optimalParams.S;
valueAtOpt = results.MinEstimatedObjective;
% optimalSampleInterval = optimalParams.D

optimalSentMat(j,k) = optimalSent/sampleSize;
optimalSentMat2(j,k) = optimalSent/maxNoSentinels;
valueAtOptMat(j,k) = valueAtOpt
end
close all
    end
end


%%
% OPTMATSAVE
% levels = [-50 -40 -30 -20 -10 -5 -1 0 5 10 20 30 40 50];
figure(1);
[C,h] = contourf(valueAtOptMat');
% [C,h] = contourf(Mstore'-M2store');


% h.LineColor = [1 0 0]; h.LineStyle = ':'; h.LineWidth = 2;
% mylabels = clabel(C,h,'manual','color','k','FontSize',16);
% for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

colbar = colorbar;
ylabel(colbar, 'Change in EDP from baseline at optimum (%)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

xlabel('Number of sentinels added (P_S)');
ylabel('Sentinel "Undetectable" period (\gamma_S)');

xticks = 1:1:length(numSentinelsVec);
xticklabels = strsplit(num2str(numSentinelsVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

yticks = 1:1:length(gammaSVec);
yticklabels = strsplit(num2str(gammaSVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

box off; set(gca,'Fontsize',16,'Linewidth',2);
%%

% CUBIC INTERPOLATION
[X,Y] = meshgrid(gammaSVec,numSentinelsVec);
gVec = 30:1:390;
sVec = 10:1:250;
[X2,Y2] = meshgrid(gVec,sVec);
Vq = interp2(X,Y,valueAtOptMat,X2,Y2,'cubic');
Vqflip = Vq';
myRowID = gammaS-29;
myRow = Vqflip(20,:);
[M,I] = min(myRow);
myRowMinID = I;
myRowMin = myRowMinID+9;
xPlot = (myRowMin+70)/80;

hold on
gammaSposition = (49+90)/120;
myline=line([1 4],[gammaSposition gammaSposition]);
myline.Color=[1 0 0]; myline.LineWidth = 2;

myMarker = plot(xPlot,gammaSposition);
myMarker.Marker = 'o'; myMarker.MarkerSize = 10;
myMarker.Color = [1 0 0]; myMarker.MarkerFaceColor = [1 0 0];

%%
% OPTMATSAVE
levels = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
figure();
[C,h] = contourf((optimalSentMat)',levels);
% [C,h] = contourf(Mstore'-M2store');


% h.LineColor = [1 0 0]; h.LineStyle = ':'; h.LineWidth = 2;
% mylabels = clabel(C,h,'manual','color','k','FontSize',16);
% for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

caxis([0 1])

colbar = colorbar;
ylabel(colbar, 'Optimal proportion of sentinels in sample', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

xlabel('Number of sentinels added (P_S)');
ylabel('Sentinel "Undetectable" period (\gamma_S)');

xticks = 1:1:length(numSentinelsVec);
xticklabels = strsplit(num2str(numSentinelsVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

yticks = 1:1:length(gammaSVec);
yticklabels = strsplit(num2str(gammaSVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

box off; set(gca,'Fontsize',16,'Linewidth',2);

% %% CUBIC INTERPOLATION
% [X,Y] = meshgrid(gammaSVec,numSentinelsVec);
% gVec = 30:1:390;
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
% % xticks = linspace(1,length(sVec),length(numSentinelsVec));
% % xticklabels = strsplit(num2str(numSentinelsVec));
% % set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);
% % 
% % yticks = linspace(1,length(gVec),length(gammaSVec));
% % yticklabels = strsplit(num2str(gammaSVec));
% % set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);
% 
% box off; set(gca,'Fontsize',16,'Linewidth',2);
% % 
% % matMin = min(min(Vq));
% % [x,y] = find(Vq==matMin);
% % optSent = sVec(x)
% % optG = gVec(y)
% 
% Vqflip = Vq';
% myRowID = gammaS-29;
% myRow = Vqflip(20,:);
% [M,I] = min(myRow);
% myRowMinID = I;
% myRowMin = myRowMinID+9;
% xPlot = (myRowMin+70)/80;

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


% leg = legend(legEntries);
% leg.Title.String = 'Asymptomatic period (days)';
% leg.FontSize = 16; leg.Location = 'northwest'; leg.Box = 'off'; 
% leg.NumColumns = 2; leg.Orientation = 'horizontal'; 


%%
function OUTPUT = testfunc(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress,baseline_EDI)
[~, EDI, ~, ~, ~, ~] = runSampling_SCI_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,"no");

OUTPUT = 100*(EDI-baseline_EDI)/baseline_EDI;
% [~, ~, ~, ~, ~, perc95] = runSampling_SCI_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,"no");
% 
% OUTPUT = 100*(perc95-baseline_95)/baseline_95;
end


