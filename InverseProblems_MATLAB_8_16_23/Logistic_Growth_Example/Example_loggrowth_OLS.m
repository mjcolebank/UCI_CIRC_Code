% In this example, we will:
% (1) Generate data from a single ODE equation and add WGN
% (2) Construct a "cost function" and perform OLS
% (3) Interpret our parameters after model calibration
clear; clc; close all;
%% A) First, we will define our equation - the Logistic Growth
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

%% B) Now, we will generate noisy, subsampled data, from our true model
data_ids = 5:5:npts;
tdata = tspace(data_ids);
ydata_clean = y(data_ids);

% Take the clean add Gaussian noise with a fixed error variance
noise_var = 5;
ydata_noisy = ydata_clean+normrnd(0,noise_var,length(data_ids),1);

% Visual the noisy data on top of the model
figure(1);clf;hold on;
plot(t,y,'k','LineWidth',3);
plot(tdata,ydata_clean,'ko','LineWidth',3,'MarkerSize',12); 
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
RSS_func = @(q) RSS_loggrowth(q,ydata_noisy,x0,tspace,data_ids);
% Now, use a built in optmization routine, like fminsearch or fmincon
par_upper = [0.8 300]; %Upper bound
par_lower = [0.01 10];  %Lower bound
par_guess = [0.1 100];   %Initial Guess
[paramopt,RSSopt] = fmincon(RSS_func,par_guess,[],[],[],[],par_lower,par_upper);
disp([param; paramopt])

% Look at the model predictions at the optimal parameter vector
f = @(t,x) logistic_growth_model(t,x,paramopt);
[t,yopt] = ode45(f,tspace,x0);

figure(1); hold on;
plot(t,yopt,'c','LineWidth',3);

% Always check your residuals!
figure(2); hold on;
plot(yopt(data_ids)-ydata_noisy,'ko','LineWidth',3,'MarkerSize',12);
grid on; set(gca,'FontSize',20);
ylabel('Residual'); xlabel('Time (days)')