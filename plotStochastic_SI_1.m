function plotStochastic_SI_1(simData,numPlots,popSize,tFinal,newFigure,yplot)

% INPUT
% simData: cell array of simulated incidence curves, generated from runSpread_SI_1.m
% numPlots: number of stochastic curves to plot
% popSize: total population size
% tFinal: final time to plot solutions for
% newFigure: specifies whether results should be plotted in a new figure ("yes") or in an
% existing figure window specified before the function call ("no")

if (newFigure ~= "yes" && newFigure ~= "no")
    fprintf('ERROR: Please enter a valid argument for newFigure ("yes" or "no")\n\n'); return
end
if (yplot ~= "num" && yplot ~= "prop")
    fprintf('ERROR: Please enter a valid argument for yplot ("num" or "prop")\n\n'); return
end

% Extract number of stochastic simulation runs available in simData
dataSize = size(simData); numSims = dataSize(1);

% If numPlots>numSims, return error message and exit
if numPlots > numSims
    fprintf(strcat('ERROR: numPlots exceeds numSims. Please set numPlots<=',num2str(numSims),'.\n\n'));
    return
end

% Choose incidence curves to plot at random
selectedRuns = randsample(numSims,numPlots);
selectedSimData = simData(selectedRuns,:);

% Plot
% Create a new figure to plot in if specified in input argument 'newFigure'
if newFigure == "yes" 
    figure(); hold on; box on; set(gca,'Fontsize',16);
end
for i = 1:numPlots
    plot(selectedSimData{i,1},selectedSimData{i,2});
end
xlim([0 tFinal]); ylim([0 popSize*1.1]);
xlabel('Time'); 
ylabel('Number of infected plants')

end