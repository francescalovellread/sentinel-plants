function [sampleData, EDI, EDT, EDR, ECDI, ESDI, perc95] = runSampling_SCI_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress)

% INPUT
% simData: cell array of simulated incidence curves, generated from runSpread_SCI_2.m
% numRuns: number of sampling runs to perform (must be <= number of simulations available in simData!)
% popSize: total population size
% numSentinels: number of sentinels in the population
% cropSampleSize: number of crop plants to sample on each sampling round
% sentinelSampleSize: number of sentinel plants to sample on each sampling round
% sampleInterval: sampling interval (time between each sampling round)
% tFinal: final permissible sample time
% progress: specifies whether progress messages are displayed ("yes" or "no")

% OUTPUT
% sampleData: a five-column matrix in which each row corresponds to a single sampling run.
% Entries in the first column are total discovery incidences; entries in the second column are
% the corresponding discovery times; entries in the third column are the numbers of sampling
% rounds that occurred. Entries in columns four and five are the crop and sentinel
% discovery incidences (i.e. total incidence broken down into species). All discovery
% incidences include both symptomatic and asymptomatic infected individuals.
% EDI: the simulated expected total discovery incidence (crops and sentinels, symptomatic
% and asymptomatic)
% EDT: the simulated expected discovery time
% EDR: the simulated expected number of sampling rounds
% ECDI: the simulated expected crop discovery incidence (symptomatic and asymptomatic)
% ESDI: the simulated expected sentinel discovery incidence (symptomatic and asymptomatic)

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

tic

Nc = cropSampleSize; Ns = sentinelSampleSize; D = sampleInterval;
P = popSize; Ps = numSentinels; Pc = P-Ps;
% Extract number of simulations from simData
dataSize = size(simData); numSims = dataSize(1);
% If numRuns>numSims, return error message and exit
if numRuns > numSims
    fprintf(strcat('Error: number of sampling runs exceeds available number of simulations. Please set numRuns<=',num2str(numSims),'.\n\n'));
    return; 
end
selectedRuns = randsample(numSims,numRuns);
selectedSimData = simData(selectedRuns,:);
% Create empty matrix to store sampling data
sampleData = zeros(numRuns,5);
% Generate initial sampling times for each run. Initial sampling times are
% uniformly distributed on the interval [0,D] to mimic pathogen
% introduction at a random time relative to the sampling scheme.
initSampleTimes = rand(numRuns,1)*D;

if progress == "yes"
    fprintf('Running sampling simulations...\t')
end
% Run sampling simulations
for i=1:numRuns
    sampleTime = initSampleTimes(i);
    detectionIndicator = 0;
    samplingRound = 1;
    while (detectionIndicator == 0 && sampleTime <= tFinal)
        sampleIndex = sum(selectedSimData{i,1} <= sampleTime); % Determine which time point on the inidence curve corresponds to the sample time
        numSymptomaticCrops = selectedSimData{i,2}(sampleIndex); % Compute total number of symptomatic crop plants
        numSymptomaticSentinels = selectedSimData{i,3}(sampleIndex); % Compute total number of symptomatic sentinel plants
        numCrypticCrops = selectedSimData{i,4}(sampleIndex); % Compute total number of cryptic/asymptomatic crop plants
        numCrypticSentinels = selectedSimData{i,5}(sampleIndex); % Compute total number of cryptic/asymptomatic sentinel plants
        cropStateVec = ones(1,numSymptomaticCrops); cropStateVec(Pc) = 0; % Create vector containing a 1 for every symptomatic crop plant and a 0 for every asymptomatic/healthy crop plant
        sentinelStateVec = ones(1,numSymptomaticSentinels); sentinelStateVec(Ps) = 0; % Create vector containing a 1 for every symptomatic sentinel plant and a 0 for every asymptomatic/healthy sentinel plant
%         selectCropVec = randsample(cropStateVec,Nc); % Take a random sample of size Nc (without replacement) from the crop vector
        selectCropVec = cropStateVec(randsample(length(cropStateVec),Nc)); % Take a random sample of size Nc (without replacement) from the crop vector
        infCropSample = sum(selectCropVec); % Total number of symptomatic crop plants in the sample
%         selectSentinelVec = randsample(sentinelStateVec,Ns); % Take a random sample of size Ns (without replacement) from the sentinel vector
        selectSentinelVec = sentinelStateVec(randsample(length(sentinelStateVec),Ns)); % Take a random sample of size Ns (without replacement) from the sentinel vector
        infSentinelSample = sum(selectSentinelVec); % Total number of symptomatic sentinel plants in the sample
             
        if infCropSample+infSentinelSample > 0 % Then an symptomatic plant has been sampled
            detectionIndicator = 1;
        else
            sampleTime = sampleTime + D; % Move on to next sample time
            samplingRound = samplingRound + 1;
        end
    end
    sampleData(i,1) = numSymptomaticCrops+numSymptomaticSentinels+numCrypticCrops+numCrypticSentinels; % Total number of infected individuals
    sampleData(i,2) = sampleTime; % Discovery time
    sampleData(i,3) = samplingRound; % Number of sampling rounds
    sampleData(i,4) = numSymptomaticCrops+numCrypticCrops; % Total number of infected crops
    sampleData(i,5) = numSymptomaticSentinels+numCrypticSentinels; % Total number of infected sentinels
end
EDI = mean(sampleData(:,1)); % Expected total discovery incidence (crops and sentinels, symptomatic and asymptomatic)
EDT = mean(sampleData(:,2)); % Expected discovery time
EDR = mean(sampleData(:,3)); % Expected number of sampling rounds
ECDI = mean(sampleData(:,4)); % Expected crop discovery incidence (symptomatic and asymptomatic)
ESDI = mean(sampleData(:,5)); % Expected sentinel discovery incidence (symptomatic and asymptomatic)
perc95 = prctile(sampleData(:,1),95);

elapsedTime = toc;
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTime),32,'secs)\n',num2str(numRuns),32,'sampling simulations performed.\n\n'));
end
end