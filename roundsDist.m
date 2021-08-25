function [roundsBars, roundsBinMidpoints] = roundsDist(sampleData,EDR,shouldIplot)

% INPUT
% sampleData: a matrix generated as the first output of a runSampling function. Each row
% corresponds to a single sampling run. Entries in the first column are discovery incidences; 
% entries in the second column are the corresponding discovery times. Entries in the third
% column are the numbers of sampling rounds before detection occurs. For two species
% models, additional fourth and fifth columns contain discovery incidences amongst crop
% and sentinel plants respectively.
% EDR: the expected number of sampling rounds from sampleData (mean of entries in third column;
% the fourth output of all runSampling functions)
% shouldIplot: a string ("yes" or "no") specifying whether the generated histogram should
% be plotted or not.

% OUTPUT
% roundsDist: vector of the normalised counts in each histogram bin
% roundsBinMidpoints: vector of the bin midpoints used to generate the histogram

if (shouldIplot ~= "yes" && shouldIplot ~= "no")
    fprintf('ERROR: Please enter a valid argument for shouldIplot ("yes" or "no")\n\n'); return
end

numSamplingRounds = sampleData(:,3); % Extract column vector of number of sampling rounds from sampleData
maxRounds = max(numSamplingRounds);

roundsBinEdges = 0.5:1:maxRounds+0.5; % Define bin edges for histogram
roundsBinMidpoints = 1:1:maxRounds; % Compute bin midpoints
roundsBars = histcounts(numSamplingRounds,roundsBinEdges); %Sort numbers of sampling rounds into bins
normaliser = sum(roundsBars); % Normalise
roundsDist = roundsBars/normaliser;

% Plot histogram of number of sampling rounds, if specified to do so in input
if shouldIplot == "yes"
    figure(); hold on; box on; set(gca,'Fontsize',16);
    roundsPlot = bar(roundsBinMidpoints,roundsDist);
    roundsPlot.FaceColor = [0.2 0.9 0.8];
    roundsPlot.EdgeColor = [0.2 0.9 0.8];
    % Plot mean number of sampling rounds
    EDRplot = xline(EDR,'k--','linewidth',3);

    xlabel('Number of sampling rounds');
    ylabel('Probability density');
    xlim([0.5 maxRounds+0.5]);
    leg = legend([roundsPlot EDRplot ],'Simulated number of sampling rounds','Expected number of sampling rounds');
    leg.Location = 'southoutside'; leg.Box = 'off'; leg.FontSize = 16;
end
end