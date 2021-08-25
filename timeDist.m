function [timeBars, timeBinMidpoints] = timeDist(sampleData,EDT,binWidth,shouldIplot)

% INPUT
% sampleData: a matrix generated as the first output of a runSampling function. Each row
% corresponds to a single sampling run. Entries in the first column are discovery incidences; 
% entries in the second column are the corresponding discovery times. Entries in the third
% column are the numbers of sampling rounds before detection occurs. For two species
% models, additional fourth and fifth columns contain discovery incidences amongst crop
% and sentinel plants respectively.
% EDT: the expected discovery time from sampleData (mean of entries in second column; the
% third ouput of all runSampling functions)
% binWidth: bin width for histogram
% shouldIplot: a string ("yes" or "no") specifying whether the generated histogram should
% be plotted or not.

% OUTPUT
% timeDist: vector of the normalised counts in each histogram bin
% timeBinMidpoints: vector of the bin midpoints used to generate the histogram

if (shouldIplot ~= "yes" && shouldIplot ~= "no")
    fprintf('ERROR: Please enter a valid argument for shouldIplot ("yes" or "no")\n\n'); return
end

discoveryTimes = sampleData(:,2); % Extract column vector of discovery times from sampleData
timeBinEdges = 0:binWidth:max(discoveryTimes)+binWidth; % Define bin edges for histogram (last+binWidth ensures last is actually contained within the final bin; otherwise the bins may end before last.
timeBinMidpoints = binWidth/2:binWidth:max(discoveryTimes)+binWidth/2; % Compute bin midpoints
timeBars = histcounts(discoveryTimes,timeBinEdges); % Sort discovery times into bins
normaliser = sum(timeBars); % Normalise
timeBars = timeBars/normaliser;

% Plot discovery time histogram, if specified to do so in input
if shouldIplot == "yes"
    figure(); hold on; box on; set(gca,'Fontsize',16);
    timesPlot = bar(timeBinMidpoints,timeBars);
    timesPlot.FaceColor = [0.7 0.9 0.6];
    timesPlot.EdgeColor = [0.7 0.9 0.6];
    % Plot mean of simulated stochastic discovery times
    EDTplot = xline(EDT,'k--','linewidth',3);

    xlabel('Discovery Time');
    ylabel('Probability density');
    xlim([0 max(discoveryTimes)]);
    leg = legend([timesPlot EDTplot],'Simulated discovery times','Expected discovery time');
    leg.Location = 'southoutside'; leg.Box = 'off'; leg.FontSize = 16;
end
end