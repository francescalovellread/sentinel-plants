function [tDet, IcDet, IsDet, ScDet, SsDet] = plotDeterministic_SI_2(popSize,numSentinels,initialI,betaCC,betaCS,betaSC,betaSS,tFinal,newFigure)

% INPUT
% popSize: total population size
% numSentinels: number of sentinels in the population
% initialI: total initial number of infected individuals (crops and sentinels)
% betaCC: transmission coefficient from crops to crops
% betaCS: transmission coefficient from crops to sentinels
% betaSC: transmission coefficient from sentinels to crops
% betaSS: transmission coefficient from sentinels to sentinels
% tFinal: final time for solution
% newFigure: specifies whether results should be plotted in a new figure ("yes") or in an
% existing figure window specified before the function call ("no")

% OUTPUT
% tDet: vector of t values for deterministic solution
% IcDet: vector of Ic (infected crop) values for deterministic solution
% IsDet: vector of Is (infected sentinel) values for deterministic solution
% ScDet: vector of Sc (susceptible crop) values for deterministic solution
% SsDet: vector of Ss (susceptible sentinel) values for deterministic solution

if (newFigure ~= "yes" && newFigure ~= "no")
    fprintf('ERROR: Please enter a valid argument for newFigure ("yes" or "no")\n\n'); return
end

P = popSize; Ps = numSentinels; Pc = P-Ps; 
bcc = betaCC; bcs = betaCS; bsc = betaSC; bss = betaSS;

% Compute proportion of sentinels in population
sentinel_prop = Ps/P; 
% Allocate the initial infected randomly between crops and sentinels according to their
% relative prevalences in the population
I0 = initialI; Ic0 = 0; Is0 = 0; rands = rand(1,I0);
for i = 1:I0
    if rands(i)<=sentinel_prop; Is0 = Is0+1;
    else 
        Ic0 = Ic0+1;
    end
end

% Compute initial numbers of susceptible crops and sentinels
Sc0 = Pc-Ic0; Ss0 = Ps-Is0; 

% Define time span
tSpan = [0 tFinal];
% Define initial conditions
f0 = [Sc0 Ss0 Ic0 Is0];
% Solve
[tDet,f] = ode45(@(tDet,f) odefcn(tDet,f,bcc,bcs,bsc,bss), tSpan, f0);
% Extract solutions
ScDet = f(:,1); SsDet = f(:,2); IcDet = f(:,3); IsDet = f(:,4);

% Plot
% Create a new figure to plot in if specified in input argument 'newFigure'
if newFigure == "yes" 
    figure();
end
hold on; box on; set(gca,'Fontsize',16);

cropIncidenceCurve = plot(tDet,IcDet);
cropIncidenceCurve.Color = [0.6 0.1 0.7]; cropIncidenceCurve.LineWidth = 2;
sentinelIncidenceCurve = plot(tDet,IsDet);
sentinelIncidenceCurve.Color = [1 0.4 0.7]; sentinelIncidenceCurve.LineWidth = 2;
xlim([0 tFinal]); ylim([0 P]);
xlabel('Time'); ylabel('Number of infected plants');
leg = legend([cropIncidenceCurve sentinelIncidenceCurve],'Infected crop plants','Infected sentinel plants');
leg.Location = 'southoutside'; leg.Box = 'off';
end

% Define ODE solver function
function dydt = odefcn(~,y,bcc,bcs,bsc,bss)
dydt = zeros(4,1);
dydt(1) = -y(1)*(bcc*y(3)+bsc*y(4));
dydt(2) = -y(2)*(bcs*y(3)+bss*y(4));
dydt(3) = y(1)*(bcc*y(3)+bsc*y(4));
dydt(4) = y(2)*(bcs*y(3)+bss*y(4));
end