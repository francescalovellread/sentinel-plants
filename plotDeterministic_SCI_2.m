function [tDet, IcDet, IsDet, CcDet, CsDet, ScDet, SsDet] = plotDeterministic_SCI_2(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,tFinal,newFigure)

% INPUT
% popSize: total population size
% numSentinels: number of sentinels in the population
% initialI: total initial number of symptomatic individuals (crops and sentinels)
% initialC: total initial number of cryptic/asymptomatic individuals (crops and sentinels)
% betaC: transmission coefficient from crops
% betaS: transmission coefficient from sentinels
% epsilonC: relative infectivity of cryptic crops compared to symptomatic crops
% epsilonS: relative infectivity of cryptic sentinels compared to symptomatic sentinels
% gammaC: cryptic/asymptomatic period for crops
% gammaS: cryptic/asymptomatic period for sentinels
% tFinal: final time for solution
% newFigure: specifies whether results should be plotted in a new figure ("yes") or in an
% existing figure window specified before the function call ("no")

% OUTPUT
% tDet: vector of t values for deterministic solution
% IcDet: vector of Ic (symptomatic crop) values for deterministic solution
% IsDet: vector of Is (symptomatic sentinel) values for deterministic solution
% CcDet: vector of Cc (cryptic crop) values for deterministic solution
% CsDet: vector of Cs (cryptic sentinel) values for deterministic solution
% ScDet: vector of Sc (susceptible crop) values for deterministic solution
% SsDet: vector of Ss (susceptible sentinel) values for deterministic solution

if (newFigure ~= "yes" && newFigure ~= "no")
    fprintf('ERROR: Please enter a valid argument for newFigure ("yes" or "no")\n\n'); return
end

P = popSize; Ps = numSentinels; Pc = P-Ps; 
bc = betaC; bs = betaS; ec = epsilonC; es = epsilonS; gc = gammaC; gs = gammaS;

% Compute proportion of sentinels in population
sentinel_prop = Ps/P; 
% Allocate the initial infected randomly between crops and sentinels according to their
% relative prevalences in the population
I0 = initialI; C0 = initialC; 
Ic0 = 0; Is0 = 0; Cc0 = 0; Cs0 = 0; 
Irands = rand(1,I0); Crands = rand(1,C0);
for i = 1:I0
    if Irands(i)<=sentinel_prop; Is0 = Is0+1;
    else 
        Ic0 = Ic0+1;
    end
end
for i = 1:C0
    if Crands(i)<=sentinel_prop; Cs0 = Cs0+1;
    else 
        Cc0 = Cc0+1;
    end
end

% Compute initial numbers of susceptible crops and sentinels
Sc0 = Pc-Ic0-Cc0; Ss0 = Ps-Is0-Cs0; 

% Define time span
tSpan = [0 tFinal];
% Define initial conditions
f0 = [Sc0 Ss0 Cc0 Cs0 Ic0 Is0];
% Solve
[tDet,f] = ode45(@(tDet,f) odefcn(tDet,f,bc,bs,ec,es,gc,gs), tSpan, f0);
% Extract solutions
ScDet = f(:,1); SsDet = f(:,2); CcDet = f(:,3); CsDet = f(:,4); IcDet = f(:,5); IsDet = f(:,6);

% Plot
% Create a new figure to plot in if specified in input argument 'newFigure'
if newFigure == "yes" 
    figure(); 
end
hold on; box on; set(gca,'Fontsize',16);

IcCurve = plot(tDet,IcDet); IcCurve.Color = [0.6 0.1 0.7]; IcCurve.LineWidth = 2;
IsCurve = plot(tDet,IsDet); IsCurve.Color = [1 0.4 0.7]; IsCurve.LineWidth = 2;
CcCurve = plot(tDet,CcDet); CcCurve.Color = [0.6 0.1 0.7]; CcCurve.LineWidth = 2; CcCurve.LineStyle = '-.';
CsCurve = plot(tDet,CsDet); CsCurve.Color = [1 0.4 0.7]; CsCurve.LineWidth = 2; CsCurve.LineStyle = '-.';
totcCurve = plot(tDet,IcDet+CcDet); totcCurve.Color = [0.6 0.1 0.7]; totcCurve.LineWidth = 0.5; totcCurve.LineStyle = ':';
totsCurve = plot(tDet,IsDet+CsDet); totsCurve.Color = [1 0.4 0.7]; totsCurve.LineWidth = 0.5; totsCurve.LineStyle = ':';

xlim([0 tFinal]); ylim([0 P]);
xlabel('Time'); ylabel('Number of infected plants');
leg = legend([IcCurve CcCurve totcCurve IsCurve CsCurve totsCurve],'Symptomatic crop plants','Asymptomatic crop plants','Total infected crop plants','Symptomatic sentinel plants','Asymptomatic sentinel plants','Total infected sentinel plants');
leg.Location = 'southoutside'; leg.Box = 'off'; leg.NumColumns = 2;
end

% Define ODE solver function
function dydt = odefcn(~,y,bc,bs,ec,es,gc,gs)
dydt = zeros(6,1);
dydt(1) = -y(1)*(bc*(ec*y(3)+y(5))+bs*(es*y(4)+y(6)));
dydt(2) = -y(2)*(bc*(ec*y(3)+y(5))+bs*(es*y(4)+y(6)));
dydt(3) = y(1)*(bc*(ec*y(3)+y(5))+bs*(es*y(4)+y(6)))-(1/gc)*y(3);
dydt(4) = y(2)*(bc*(ec*y(3)+y(5))+bs*(es*y(4)+y(6)))-(1/gs)*y(4);
dydt(5) = (1/gc)*y(3);
dydt(6) = (1/gs)*y(4);
end