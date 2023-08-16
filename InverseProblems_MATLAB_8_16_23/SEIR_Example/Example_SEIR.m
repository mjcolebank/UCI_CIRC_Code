% In this example, we will look at a multistate model describing epidemic
% dynamics
clear; clc; close all;
%% First, we will define our equation - the SEIR model
%       dS/dt = -beta.*S.*I./N - mu.*S + mu.*N; (susceptible)
%       dE/dt = beta.*S.*I./N - (mu+delta(t)).*E;  (exposed)
%       dI/dt = delta(t).*E - (gamma+mu).*I;       (infected)
%       dR/dt = gamma.*I - mu.*R;               (recovered)

% NOTE: delta represents the exposed->infected rate. We can make this a
% 'seasonal' parameter by setting 'infect_cycle>0'. Otherwise, setting
% infect_cycle = 0 gives us delta(t) = fixed.

% Define the abscissa for the model (here, time in days)
npts   = 101;
Tend = 150;
tspace = linspace(0,Tend,npts); % This will be our time interval in days

% Define our initial conditions  for the ODE model and the parameters
S0 = 999; E0 = 0; I0 = 1; R0 = 0;
beta         = 0.8; % infection rate
mu           = 0.1; % birth/death rate
gamma        = 0.1; % recovery after infection rate
delta        = 0.4; % exposed to infected rate (could be seasonal)
infect_cycle = 25;  % describes when infection occurs in time

% To make things easier, stack our parameter values into a vector
param = [beta mu gamma delta infect_cycle];
x0 = [S0; E0; I0; R0]; N = sum(x0);
% Visualize the model output and its derivative
f = @(t,x) SEIR_model(t,x,param,N);
[t,y] = ode45(f,tspace,x0);
dy = zeros(npts,4);
for i=1:npts
    dy(i,:)    = f(t(i),y(i,:)); %Evaluate the right handside for plotting purposes
end
%%
% Plot the model and the derivative (rhs of ODE)
figure(1);clf;
plot(t,y,'LineWidth',3); ylabel('States'); 
grid on; set(gca,'FontSize',20); axis tight; xlabel('Time (days)');
legend('S','E','I','R')

figure(2);clf;
plot(t,dy,'LineWidth',3); ylabel('dy/dt'); 
grid on; set(gca,'FontSize',20); axis tight; xlabel('Time (days)');
