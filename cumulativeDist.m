function [cdfBars, cdfBinMidpoints] = cumulativeDist(sampleData,popSize,plotStyle)

% INPUT
% sampleData: a matrix generated as the first output of a runSampling function. Each row
% corresponds to a single sampling run. Entries in the first column are discovery incidences; 
% entries in the second column are the corresponding discovery times. Entries in the third
% column are the numbers of sampling rounds before detection occurs. For two species
% models, additional fourth and fifth columns contain discovery incidences amongst crop
% and sentinel plants respectively.
% popSize: the total population size
% plotStyle: a string ("bars" or "line") specifying whether the resulting CDF should be
% plotted as a bar chart or a line graph

% OUTPUT
% cumulativeDist: vector of the normalised CDF counts in each histogram bin
% binMidpoints: vector of the bin midpoints used to generate the histogram

if (plotStyle ~= "bars" && plotStyle ~= "line")
    fprintf('ERROR: Please enter a valid argument for plotStyle ("bars" or "line")\n\n'); return
end
        
discoveryIncidences = sampleData(:,1); % Extract column vector of discovery incidences from sampleData
binEdges = 0.5:1:popSize+0.5; % Define bin edges for histogram
cdfBinMidpoints = 1:1:popSize; % Compute bin midpoints
incidenceBars = histcounts(discoveryIncidences,binEdges); % Sort discovery incidences into bins
normaliser = sum(incidenceBars); % Normalise
incidenceDist = incidenceBars/normaliser;
cdfBars = cumsum(incidenceDist);

% Plot discovery incidence histogram, if specified to do so in input
figure(); hold on; box on; set(gca,'Fontsize',16);
if  plotStyle == "bars"
    cumulativePlot = bar(cdfBinMidpoints, cdfBars);
    cumulativePlot.FaceColor = [0.2 0.7 1];
    cumulativePlot.EdgeColor = [0.2 0.7 1];
elseif plotStyle == "line"
    cumulativePlot = plot(cdfBinMidpoints, cdfBars);
    cumulativePlot.Color = [0.2 0.7 1];
    cumulativePlot.LineWidth = 2;
end
limline = line([0,popSize],[1,1]); limline.LineStyle='--';
xlabel('Discovery incidence (I^*)');
ylabel('Cumulative distribution');
xlim([0, popSize+0.5]); ylim([0 1.1]);
end