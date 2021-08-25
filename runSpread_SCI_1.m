function simData = runSpread_SCI_1(popSize,initialI,initialC,beta,epsilon,gamma,numSims,tFinal,progress)

% INPUT
% popSize: total population size.
% initialI: initial number of symptomatic infected individuals
% initialC: initial number of cryptic/asymptomatic infected individuals
% beta: transmission coefficient for symptomatic individuals
% epsilon: transmission coefficient scaling factor for cryptic/asymptomatic individuals
% (e.g. if cryptic hosts are half as infectious as symptomatic hosts, set epsilon=0.5)
% asympPeriod: duration of asymptomatic period
% numSims: number of simulations to run
% tFinal: maximum time for each simulation to run
% progress: specifies whether progress messages are displayed ("yes" or "no")

% OUTPUT
% sim_data: a (no_sims x 4) cell array containing simulation results. Each
% row corresponds to a run; first, second, third and fourth columns contain
% vectors for t, I, C and S respectively.

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

tic

P = popSize; I0 = initialI; C0 = initialC; b = beta; e = epsilon; g = gamma;
% Compute initial number of susceptible individuals
S0 = P-I0-C0;
% Create cell array for storing results
simData = cell(numSims,4);
if progress == "yes"
    fprintf('Running spread simulations...\t')
end
% Run stochatic simulations
for i=1:numSims
    t = 0; I = I0; C = C0; S = S0;
    tvec = zeros(1,1+2*S0+C0); Ivec = zeros(1,1+2*S0+C0); Cvec = zeros(1,1+2*S0+C0); Svec = zeros(1,1+2*S0+C0);
    index = 1; tvec(index)=t; Ivec(index)=I0; Cvec(index)=C0; Svec(index)=S0;
    while (t<tFinal && (S+C)>0)
        a1 = b*S*I; a2 = b*e*S*C; a3 = (1/g)*C; % Compute individual reaction propensities
        a0 = a1+a2+a3; % Compute total reaction propensity
        r1 = rand(1); tau = (1/a0)*log(1/r1); % Compute time to next reaction
        r2 = rand(1);
        if r2<(a1+a2)/a0 % Then an S becomes a C
            C = C+1; S = S-1;
        else % Then a C becomes an I
            I = I+1; C = C-1;
        end
        t = t+tau;
        index = index + 1; 
        tvec(index) = t; 
        Ivec(index) = I; Cvec(index) = C; Svec(index) = S; 
    end
    simData{i,1} = tvec;
    simData{i,2} = Ivec;
    simData{i,3} = Cvec;
    simData{i,4} = Svec;
end
elapsedTime = toc;
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTime),32,'secs)\n',num2str(numSims),32,'incidence curves generated.\n\n'));
end
end
