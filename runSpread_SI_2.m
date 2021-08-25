function simData = runSpread_SI_2(popSize,numSentinels,initialI,betaCC,betaCS,betaSC,betaSS,numSims,tFinal,progress)

% INPUT
% popSize: total population size.
% numSentinels: number of sentinels in the population
% initialI: total initial number of infected individuals (crops and sentinels)
% betaCC: transmission coefficient from crops to crops
% betaCS: transmission coefficient from crops to sentinels
% betaSC: transmission coefficient from sentinels to crops
% betaSS: transmission coefficient from sentinels to sentinels
% numSims: number of simulations to run
% tFinal: maximum time for each simulation to run
% progress: specifies whether progress messages are displayed ("yes" or "no")

% OUTPUT
% sim_data: a (no_sims x 5) cell array containing simulation results. Each
% row corresponds to a run; first, second, third fourth and fifth columns
% contain vectors for t, Ic (infected crops), Is (infected sentinels), Sc
% (susceptible crops) and Ss (susceptipble sentinels) respectively.

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

tic

P = popSize; Ps = numSentinels; Pc = P-Ps; sentinel_prop = Ps/P;
bcc = betaCC; bcs = betaCS; bsc = betaSC; bss = betaSS; 
% Create cell array for storing results
simData = cell(numSims,5);
if progress == "yes"
    fprintf('Running spread simulations...\t')
end

% Run stochatic simulations
for i=1:numSims
    % On each run, divide initial infected individuals between crops and sentinels
    % according to their relative proportions in the population
    I0 = initialI; Ic0 = 0; Is0 = 0; rands = rand(1,I0);
    for j = 1:I0
        if rands(j)<=sentinel_prop
            Is0 = Is0+1;
        else
            Ic0 = Ic0+1;
        end
    end
    % Compute initial numbers of susceptible crops and sentinels
    Sc0 = Pc-Ic0; Ss0 = Ps-Is0;
    
    t = 0; Ic = Ic0; Is = Is0; Sc = Sc0; Ss = Ss0;
    tvec = zeros(1,1+Sc0+Ss0); Icvec = zeros(1,1+Sc0+Ss0); Isvec = zeros(1,1+Sc0+Ss0); Scvec = zeros(1,1+Sc0+Ss0); Ssvec = zeros(1,1+Sc0+Ss0);
    index = 1; tvec(index)=t; Icvec(index)=Ic0; Isvec(index)=Is0; Scvec(index)=Sc0; Ssvec(index)=Ss0;
    while (t<tFinal && Sc+Ss>0)
        a1 = bcc*Sc*Ic; a2 = bsc*Sc*Is; a3 = bcs*Ss*Ic; a4 = bss*Ss*Is; % Compute individual reaction propensities
        a0 = a1+a2+a3+a4; % Compute total reaction propensity
        r1 = rand(1); tau = (1/a0)*log(1/r1); % Computd time to next reaction
        r2 = rand(1);
        if r2<=(a1+a2)/a0 % Then a susceptible crop becomes an infected crop
            Sc = Sc-1; Ic = Ic+1;
        else % A susceptible sentinel becomes an infected sentinel
            Ss = Ss-1; Is = Is+1;
        end
        t = t+tau;
        index = index + 1;
        tvec(index) = t;
        Icvec(index) = Ic; Isvec(index) = Is; Scvec(index) = Sc; Ssvec(index) = Ss;
    end
    simData{i,1} = tvec;
    simData{i,2} = Icvec;
    simData{i,3} = Isvec;
    simData{i,4} = Scvec;
    simData{i,5} = Ssvec;
end
elapsedTime = toc;
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTime),32,'secs)\n',num2str(numSims),32,'incidence curves generated.\n\n'));
end
end
