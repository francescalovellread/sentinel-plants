function [totalIncidenceBars, incBinMidpoints] = incidenceDistMultipleSpecies(sampleData,EDI,ECDI,ESDI,first,last,binWidth,shouldIplot)

% INPUT
% sampleData: a matrix generated as the first output of a runSampling function. Each row
% corresponds to a single sampling run. Entries in the first column are discovery incidences; 
% entries in the second column are the corresponding discovery times. Entries in the third
% column are the numbers of sampling rounds before detection occurs. For two species
% models, additional fourth and fifth columns contain discovery incidences amongst crop
% and sentinel plants respectively.
% EDI: the expected discovery incidence from sampleData (mean of entries in first column; 
% the second output of all runSampling functions)
% ECDI: the expected crop discovery incidence from sampleData (mean of entries in fourth column; 
% the fifth output of two species runSampling functions)
% ESDI: the expected sentinel discovery incidence from sampleData (mean of entries in fifth column; 
% the sixth output of two species runSampling functions)
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

totalDiscoveryIncidences = sampleData(:,1); % Extract column vector of discovery incidences from sampleData
cropDiscoveryIncidences = sampleData(:,4); % Extract column vector of crop discovery incidences from sampleData
sentinelDiscoveryIncidences = sampleData(:,5); % Extract column vector of sentinel discovery incidences from sampleData

binEdges = first:binWidth:last+binWidth; % Define bin edges for histogram (last+binWidth ensures last is actually contained within the final bin; otherwise the bins may end before last.
incBinMidpoints = first+binWidth/2:binWidth:last+binWidth/2; % Compute bin midpoints

totalIncidenceBars = histcounts(totalDiscoveryIncidences,binEdges); % Sort discovery incidences into bins
normaliser = sum(totalIncidenceBars); % Normalise
totalIncidenceBars = totalIncidenceBars/normaliser;

cropIncidenceBars = histcounts(cropDiscoveryIncidences,binEdges); % Sort discovery incidences into bins
normaliser = sum(cropIncidenceBars); % Normalise
cropIncidenceBars = cropIncidenceBars/normaliser;

sentinelIncidenceBars = histcounts(sentinelDiscoveryIncidences,binEdges); % Sort discovery incidences into bins
normaliser = sum(sentinelIncidenceBars); % Normalise
sentinelIncidenceBars = sentinelIncidenceBars/normaliser;

% Plot discovery incidence histogram, if specified to do so in input
if shouldIplot == "yes"
    figure(); hold on; box on; set(gca,'Fontsize',16);
    totalIncidencePlot = bar(incBinMidpoints,totalIncidenceBars);
    totalIncidencePlot.FaceColor = [0.2 0.7 1];
    totalIncidencePlot.EdgeColor = [0.2 0.7 1];
    % Plot mean of simulated stochastic discovery incidences
    EDIplot = xline(EDI,'k--','linewidth',3);

    xlabel('Discovery incidence (I^*)');
    ylabel('Probability density');
    xlim([0, last]);
    leg = legend([totalIncidencePlot EDIplot],'Simulated discovery incidences','Expected discovery incidence');
    leg.Location = 'southoutside'; leg.Box = 'off'; leg.FontSize = 16;

    figure(); hold on; box on; set(gca,'Fontsize',16);
    cropIncidencePlot = bar(incBinMidpoints,cropIncidenceBars);
    cropIncidencePlot.FaceColor = [0.2 0.7 1];
    cropIncidencePlot.EdgeColor = [0.2 0.7 1];
    % Plot mean of simulated stochastic crop discovery incidences
    ECDIplot = xline(ECDI,'k--','linewidth',3);

    xlabel('Crop discovery incidence');
    ylabel('Probability density');
    xlim([0, last]);
    leg = legend([cropIncidencePlot ECDIplot],'Simulated crop discovery incidences','Expected crop discovery incidence');
    leg.Location = 'southoutside'; leg.Box = 'off'; leg.FontSize = 16;

    figure(); hold on; box on; set(gca,'Fontsize',16);
    sentinelIncidencePlot = bar(incBinMidpoints,sentinelIncidenceBars);
    sentinelIncidencePlot.FaceColor = [0.2 0.7 1];
    sentinelIncidencePlot.EdgeColor = [0.2 0.7 1];
    % Plot mean of simulated stochastic discovery incidences
    ESDIplot = xline(ESDI,'k--','linewidth',3);

    xlabel('Sentinel discovery incidence');
    ylabel('Probability density');
    xlim([0, last]);
    leg = legend([sentinelIncidencePlot ESDIplot],'Simulated sentinel discovery incidences','Expected sentinel discovery incidence');
    leg.Location = 'southoutside'; leg.Box = 'off'; leg.FontSize = 16;
end
end