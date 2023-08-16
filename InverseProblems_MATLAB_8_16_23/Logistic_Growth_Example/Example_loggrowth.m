% In this example, we will look at a simple ODE model of logistic growth
% that can often be used to describe ecological or biological organism
% growth
clear; clc; close all;
%% First, we will define our equation - the Logistic Growth
%       dx/dt = r*x*(1-x/K)
% where r (1/days) is the growth rate and K (# of organisms) is the carrying
% capacity

% Define the abscissa for the model (here, time in days)
npts   = 101;
tspace = linspace(0,100,npts); % This will be our time interval in days

% Define our initial conditions  for the ODE model and the parameters
x0 = 2;   % Initial amount of organisms (e.g., cells)
r = 0.12; % Rate of growth 
K = 150;  % Carrying capacity due to space, resources, etc

% To make things easier, stack our parameter values into a vector
param = [r K];

% Visualize the model output and its derivative
f = @(t,x) logistic_growth_model(t,x,param);
[t,y] = ode45(f,tspace,x0);
dy    = f(tspace,y); %Evaluate the right handside for plotting purposes

% Plot the model and the derivative (rhs of ODE)
figure(1);clf;
plot(t,y,'LineWidth',3); ylabel('Organism'); axis tight
yyaxis('right'); 
plot(t,dy,'LineWidth',3); ylabel('Rate of change');
grid on; set(gca,'FontSize',20); axis tight; xlabel('Time (days)')

