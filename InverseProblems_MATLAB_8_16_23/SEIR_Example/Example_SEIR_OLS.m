% In this example, we will:
% (1) Generate data from a single ODE equation and add WGN
% (2) Construct a "cost function" and perform OLS
% (3) Interpret our parameters after model calibration
clear; clc; close all;
%% A) First, we will define our equations
npts   = 101;
Tend = 150;
tspace = linspace(0,Tend,npts); % This will be our time interval in days

% Define our initial conditions  for the ODE model and the parameters
S0 = 999; E0 = 0; I0 = 1; R0 = 0;
beta         = 0.8; % 
mu           = 0.1; %
gamma        = 0.1; %
delta        = 0.4; %
infect_cycle = 25; % 

% To make things easier, stack our parameter values into a vector
param = [beta mu gamma delta infect_cycle];
x0 = [S0; E0; I0; R0]; N = sum(x0);

% Visualize the model output and its derivative
f = @(t,x) SEIR_model(t,x,param,N);
[t,y] = ode45(f,tspace,x0);

%% B) Now, we will generate noisy, subsampled data, from our true model
data_ids = 5:10:npts;
tdata = tspace(data_ids);
ydata_clean = y(data_ids,3);

% Take the clean add Gaussian noise with a fixed error variance
noise_var = 5;
noise_mean = 10;
ydata_noisy = ydata_clean+normrnd(noise_mean,noise_var,length(data_ids),1);
% Make sure signal is positive
ydata_noisy = max(ydata_noisy,0);
% Visual the noisy data on top of the model
figure(1);clf;hold on;
plot(t,y(:,3),'k','LineWidth',1);
plot(tdata,ydata_clean,'ko','LineWidth',3,'MarkerSize',8); 
plot(tdata,ydata_noisy,'--r+','LineWidth',3,'MarkerSize',12); 
ylabel('Organism'); xlabel('Time (days)');
grid on; set(gca,'FontSize',20);  axis tight; 

%% C) Consider model calibration
% In OLS, we want to minimize a cost-function OR maximize a likelihood
% (more on this later)
% Here, lets define the residual sum of squares
%       J = sum( (data-model).^2)
% or, in vector notation
%       J = (data-model)'*(data-model)
% This is easier to do in a separate function (see RSS_loggrowth.m)

% Create an in-line function handle 
% NOTE: we don't estimate the infection time here, so just add it as the
% last component of the vector 'q'
RSS_func = @(q) RSS_SEIR([q infect_cycle],ydata_noisy,x0,tspace,data_ids);

% Now, use a built in optmization routine, like fminsearch or fmincon
par_upper = [0.95 0.95 0.95 0.95]; %Upper bound
par_lower = [0.01 0.01 0.01 0.01];  %Lower bound
par_guess = [0.75 0.05 0.15 0.3];   %Initial Guess
[paramopt,RSSopt] = fmincon(RSS_func,par_guess,[],[],[],[],par_lower,par_upper);
disp([param(1:4); paramopt]')
%%
% Look at the model predictions at the optimal parameter vector
f = @(t,x) SEIR_model(t,x,[paramopt infect_cycle],N);
[t,yopt] = ode45(f,tspace,x0);

figure(1); hold on;
plot(t,yopt(:,3),'c','LineWidth',3);

% Always check your residuals!
figure(2); hold on;
plot(yopt(data_ids,3)-ydata_noisy,'ko','LineWidth',3,'MarkerSize',12);
grid on; set(gca,'FontSize',20);
ylabel('Residual'); xlabel('Time (days)')