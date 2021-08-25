function [tDet, IDet, CDet, SDet] = plotDeterministic_SCI_1(popSize,initialI,initialC,beta,epsilon,gamma,tFinal,newFigure,yplot)

% INPUT
% popSize: total population size
% initialI: initial number of symptomatic individuals
% initialC: initial number of cryptic/asymptomatic individuals
% betaI: transmission coefficient for infections coming from symptomatic individuals
% betaC: transmission coefficient for infections coming from cryptic/asymptomatic individuals
% gamma: duration of cryptic/asymptomatic period
% tFinal: final time for solution
% newFigure: specifies whether results should be plotted in a new figure ("yes") or in an
% existing figure window specified before the function call ("no")

% OUTPUT
% tDet: vector of t values for deterministic solution
% IDet: vector of I values for deterministic solution
% CDet: vector of C values for deterministic solution
% SDet: vector of S values for deterministic solution

if (newFigure ~= "yes" && newFigure ~= "no")
    fprintf('ERROR: Please enter a valid argument for newFigure ("yes" or "no")\n\n'); return
end
if (yplot ~= "num" && yplot ~= "prop")
    fprintf('ERROR: Please enter a valid argument for yplot ("num" or "prop")\n\n'); return
end

P = popSize; I0 = initialI; C0 = initialC; S0 = P-I0-C0; b = beta; e = epsilon; g = gamma;

% Define time span
tSpan = [0 tFinal];
% Define initial conditions
f0 = [S0 C0 I0];
% Solve
[tDet,f] = ode45(@(tDet,f) odefcn(tDet,f,b,e,g), tSpan, f0);
% Extract solutions
SDet = f(:,1); CDet = f(:,2); IDet = f(:,3);

% Plot
% Create a new figure to plot in if specified in input argument 'newFigure'
if newFigure == "yes" 
    figure(); hold on; box on; set(gca,'Fontsize',16);
end
if yplot == "num"
    incidenceCurveI = plot(tDet,IDet);
    incidenceCurveC = plot(tDet,CDet);
    incidenceCurveTot = plot(tDet,IDet+CDet);
    ylim([0 P]); ylabel('Number of infected plants');
else
    incidenceCurveI = plot(tDet,IDet/popSize);
    incidenceCurveC = plot(tDet,CDet/popSize);
    incidenceCurveTot = plot(tDet,(IDet+CDet)/popSize);
    ylim([0 1]); ylabel('Proportion of infected plants');
end

incidenceCurveI.Color = [0.6 0.1 0.7]; incidenceCurveI.LineWidth = 2;
incidenceCurveC.Color = [0.5 0.9 0.7]; incidenceCurveC.LineWidth = 2;
incidenceCurveTot.Color = 'k'; incidenceCurveTot.LineStyle = '--';

xlim([0 tFinal]); 
xlabel('Time'); 

leg = legend([incidenceCurveI incidenceCurveC incidenceCurveTot],'Symptomatic infected plants','Asymptomatic infected plants','Total infected plants');
leg.Location = 'southoutside'; leg.Box = 'off';

end

% Define ODE solver function
function dydt = odefcn(~,y,b,e,g)
dydt = zeros(3,1);
dydt(1) = -(b*y(3)+b*e*y(2))*y(1);
dydt(2) = (b*y(3)+b*e*y(2))*y(1) - (1/g)*y(2);
dydt(3) = (1/g)*y(2);
end