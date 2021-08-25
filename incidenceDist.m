function [incidenceBars, incBinMidpoints] = incidenceDist(sampleData,EDI,first,last,binWidth,shouldIplot)

% INPUT
% sampleData: a matrix generated as the first output of a runSampling function. Each row
% corresponds to a single sampling run. Entries in the first column are discovery incidences; 
% entries in the second column are the corresponding discovery times. Entries in the third
% column are the numbers of sampling rounds before detection occurs. For two species
% models, additional fourth and fifth columns contain discovery incidences amongst crop
% and sentinel plants respectively.
% EDI: the expected discovery incidence from sampleData (mean of entries in first column; 
% the second output of all runSampling functions)
% first: lower limit for histogram bin edges
% last: upper limit for histogram bin edges
% binWidth: bin width for histogram
% shouldIplot: a string ("yes" or "no") specifying whether the generated histogram should
% be plotted or not.

% OUTPUT
% incidenceDist: vector of the normalised counts in each histogram bin
% binMidpoints: vector of the bin midpoints used to generate the histogram

if (shouldIplot ~= "yes" && shouldIplot ~= "no")
    fprintf('ERROR: Please enter a valid argument for shouldIplot ("yes" or "no")\n\n'); return
end

discoveryIncidences = sampleData(:,1); % Extract column vector of discovery incidences from sampleData
binEdges = first:binWidth:last+binWidth; % Define bin edges for histogram (last+binWidth ensures last is actually contained within the final bin; otherwise the bins may end before last.
incBinMidpoints = first+binWidth/2:binWidth:last+binWidth/2; % Compute bin midpoints
incidenceBars = histcounts(discoveryIncidences,binEdges); % Sort discovery incidences into bins
normaliser = sum(incidenceBars); % Normalise
incidenceBars = incidenceBars/normaliser;

% Plot discovery incidence histogram, if specified to do so in input
if shouldIplot == "yes"
    figure(); hold on; box on; set(gca,'Fontsize',16);
    incidencePlot = bar(incBinMidpoints,incidenceBars);
    incidencePlot.FaceColor = [0.2 0.7 1];
    incidencePlot.EdgeColor = [0.2 0.7 1];
    % Plot mean of simulated stochastic discovery incidences
    EDIplot = xline(EDI,'k--','linewidth',3);

    xlabel('Discovery incidence (I^*)');
    ylabel('Probability density');
    xlim([first, last]);
    leg = legend([incidencePlot EDIplot],'Simulated discovery incidences','Expected discovery incidence');
    leg.Location = 'southoutside'; leg.Box = 'off'; leg.FontSize = 16;
end
end