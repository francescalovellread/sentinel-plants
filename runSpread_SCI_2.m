function simData = runSpread_SCI_2(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress)

% INPUT
% popSize: total population size.
% numSentinels: number of sentinels in the population
% initialI: total initial number of symptomatic individuals (crops and sentinels)
% initialC: total inital number of cryptic/asymptomatic individuals (crops and sentinels)
% betaC: transmission coefficient for symptomatic crops
% betaS: transmission coefficient for symptomatic sentinels
% epsilonC: transmission coefficient scaling factor for cryptic/asymptomatic crops
% epsilonS: transmission coefficient scaling factor for cryptic/asymptomatic sentinels
% gammaC: cryptic/asymptomatic period for crops
% gammaS: cryptic/asypmtomatic period for sentinels
% numSims: number of simulations to run
% tFinal: maximum time for each simulation to run
% progress: specifies whether progress messages are displayed ("yes" or "no")

% OUTPUT
% sim_data: a (no_sims x 7) cell array containing simulation results. Each row corresponds
% to a run; columns contain vectors for t, Ic (symptomatic crops), Is (symptomatic sentinels),
% Cc (cryptic crops), Cs (cryptic sentinels), Sc (susceptible crops) and Ss (susceptible
% sentinels).

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

tic

P = popSize; Ps = numSentinels; Pc = P-Ps; sentinel_prop = Ps/P;
bc = betaC; bs = betaS; ec = epsilonC; es = epsilonS; gc = gammaC; gs = gammaS;
% Create cell array for storing results
simData = cell(numSims,7);
if progress == "yes"
    fprintf('Running spread simulations...\t')
end

% Run stochatic simulations
for i=1:numSims
    % On each run, divide initial infected individuals between crops and sentinels
    % according to their relative proportions in the population
    I0 = initialI; C0 = initialC;
    Ic0 = 0; Is0 = 0; Cc0 = 0; Cs0 = 0; randsI = rand(1,I0); randsC = rand(1,C0);
    for j = 1:I0
        if randsI(j)<=sentinel_prop
            Is0 = Is0+1;
        else
            Ic0 = Ic0+1;
        end
    end
    for j = 1:C0
        if randsC(j)<=sentinel_prop
            Cs0 = Cs0+1;
        else
            Cc0 = Cc0+1;
        end
    end
    % Compute initial numbers of susceptible crops and sentinels
    Sc0 = Pc-Ic0-Cc0; Ss0 = Ps-Is0-Cs0;
    
    t = 0; Ic = Ic0; Is = Is0; Cc = Cc0; Cs = Cs0; Sc = Sc0; Ss = Ss0;
    vecLength = 1+2*(Sc0+Ss0)+Cc0+Cs0;
    tvec = zeros(1,vecLength);
    Icvec = zeros(1,vecLength); Isvec = zeros(1,vecLength);
    Ccvec = zeros(1,vecLength); Csvec = zeros(1,vecLength);
    Scvec = zeros(1,vecLength); Ssvec = zeros(1,vecLength);
    index = 1; tvec(index)=t; Icvec(index)=Ic0; Isvec(index)=Is0; Ccvec(index)=Cc0; Csvec(index)=Cs0; Scvec(index)=Sc0; Ssvec(index)=Ss0;
   
    while (t<tFinal && Sc+Ss+Cc+Cs>0)
        % Compute individual reaction propensities
        a1 = Sc*(bc*Ic+bs*Is+ec*bc*Cc+es*bs*Cs);
        a2 = Ss*(bc*Ic+bs*Is+ec*bc*Cc+es*bs*Cs);
        a3 = (1/gc)*Cc;
        a4 = (1/gs)*Cs;
        a0 = a1+a2+a3+a4; % Compute total reaction propensity
        r1 = rand(1); tau = (1/a0)*log(1/r1); % Compute time to next reaction
        r2 = rand(1);
        if r2<=a1/a0 % Then a susceptible crop becomes a cryptic crop
            Sc = Sc-1; Cc = Cc+1;
        elseif r2<=(a1+a2)/a0 % Then a susceptible sentinel becomes a cryptic sentinel
            Ss = Ss-1; Cs = Cs+1;
        elseif r2<=(a1+a2+a3)/a0 % Then a cryptic crop becomes a symptomatic crop
            Cc = Cc-1; Ic = Ic+1;
        else % A cryptic sentinel becomes a symptomatic sentinel
            Cs = Cs-1; Is = Is+1;
        end
        t = t+tau;
        index = index + 1;
        tvec(index) = t;
        Icvec(index) = Ic; Isvec(index) = Is;
        Scvec(index) = Sc; Ssvec(index) = Ss;
        Ccvec(index) = Cc; Csvec(index) = Cs;
    end
    simData{i,1} = tvec;
    simData{i,2} = Icvec; simData{i,3} = Isvec;
    simData{i,4} = Ccvec; simData{i,5} = Csvec;
    simData{i,6} = Scvec; simData{i,7} = Ssvec;
end
elapsedTime = toc;
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTime),32,'secs)\n',num2str(numSims),32,'incidence curves generated.\n\n'));
end
end
