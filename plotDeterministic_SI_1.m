function [tDet, IDet, SDet] = plotDeterministic_SI_1(popSize,initialI,beta,tFinal,newFigure,yplot)

% INPUT
% popSize: total population size
% initialI: initial number of infected individuals
% beta: transmission coefficient
% tFinal: final time for solution
% newFigure: specifies whether results should be plotted in a new figure ("yes") or in an
% existing figure window specified before the function call ("no")

% OUTPUT
% tDet: vector of t values for deterministic solution
% IDet: vector of I values for deterministic solution
% SDet: vector of S values for deterministic solution

if (newFigure ~= "yes" && newFigure ~= "no")
    fprintf('ERROR: Please enter a valid argument for newFigure ("yes" or "no")\n\n'); return
end
if (yplot ~= "num" && yplot ~= "prop")
    fprintf('ERROR: Please enter a valid argument for yplot ("num" or "prop")\n\n'); return
end

P = popSize; I0 = initialI; S0 = P-I0; b = beta;

% Define time span
tSpan = [0 tFinal];
% Define initial conditions
f0 = [S0 I0];
% Solve
[tDet,f] = ode45(@(tDet,f) odefcn(tDet,f,b), tSpan, f0);
% Extract solutions
SDet = f(:,1); IDet = f(:,2);

% Plot
% Create a new figure to plot in if specified in input argument 'newFigure'
if newFigure == "yes" 
    figure(); hold on; box on; set(gca,'Fontsize',16);
end
if yplot == "num"
    incidenceCurve = plot(tDet,IDet); ylim([0 P]); ylabel('Number of infected plants');
else incidenceCurve = plot(tDet,IDet/popSize); ylim([0 1]); ylabel('Proportion of infected plants');
end
incidenceCurve.Color = [0.6 0.1 0.7]; incidenceCurve.LineWidth = 2;
xlim([0 tFinal]); 
xlabel('Time'); 

leg = legend([incidenceCurve],'Total infected plants');
leg.Location = 'southoutside'; leg.Box = 'off';

end

% Define ODE solver function
function dydt = odefcn(~,y,b)
dydt = zeros(2,1);
dydt(1) = -b*y(1)*y(2);
dydt(2) = b*y(1)*y(2);
end