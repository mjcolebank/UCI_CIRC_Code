% Try running Metropolis Hastings with the Arterial wall model
clear; clc; close all;
addpath('../Algorithms/')
param_Names = {'\beta','\mu','\gamma','\delta'};
%%
npts   = 101;
Tend = 150;
tspace = linspace(0,Tend,npts); % This will be our time interval in days

% Define our initial conditions  for the ODE model and the parameters
S0 = 999; E0 = 0; I0 = 1; R0 = 0;
beta = 0.8; mu = 0.1; gamma = 0.1; delta = 0.4; infect_cycle = 25; % 

% To make things easier, stack our parameter values into a vector
param_true = [beta mu gamma delta infect_cycle];
x0 = [S0; E0; I0; R0]; N = sum(x0);

% Visualize the model output and its derivative
f = @(t,x) SEIR_model(t,x,param_true,N);
[t,ytrue] = ode45(f,tspace,x0);

%% B) Now, we will generate noisy, subsampled data, from our true model
data_ids = 5:5:npts; tdata = tspace(data_ids); ydata_clean = ytrue(data_ids,3);

% Take the clean add Gaussian noise with a fixed error variance
noise_var = 6; noise_mean = 8;
ydata_noisy = ydata_clean+normrnd(noise_mean,noise_var,length(data_ids),1);
% Make sure signal is positive
ydata_noisy = max(ydata_noisy,0);


%% 
par_upper = [0.95 0.95 0.95 0.95]; %Upper bound
par_lower = [0.01 0.01 0.01 0.01];  %Lower bound
par_guess = [0.75 0.05 0.15 0.3];   %Initial Guess
RSS_func = @(q) RSS_SEIR([q infect_cycle],ydata_noisy,x0,tspace,data_ids);

% If possible, optimize first
[par_guess,RSSopt] = fmincon(RSS_func,par_guess,[],[],[],[],par_lower,par_upper);

N_par = length(par_guess);
N_data = length(ydata_noisy);
%%
burnin = 1000;
noise_est = RSS_SEIR([par_guess infect_cycle],ydata_noisy,x0,tspace,data_ids)./(N_data-N_par);
% noise_est = 500;
inputs.par0 = par_guess;
inputs.nsamp = 10000;
inputs.ss_fun = @(q,data) RSS_SEIR([q infect_cycle],data,x0,tspace,data_ids);
inputs.data   = ydata_noisy;
inputs.measurement_noise = noise_est;
inputs.par_low = par_lower;
inputs.par_upp = par_upper;
inputs.burnin = burnin;

% If you want a Gaussian prior, uncomment
% inputs.prior = 'Gaussian';
% inputs.thetamu  = par_guess;
% inputs.thetasig = 0.5.*par_guess;

% Choose between Metropolis-Hastings or Adaptive Metropolis
chain = MH_Algorithm(inputs);
% chain = AM_Algorithm(inputs);

% Throw away the "burnin" period
chain = chain(burnin+1:end,:);

%% Chains
figure(1);clf;
for i=1:N_par
subplot(2,2,i); hold on;
plot(chain(:,i));
plot(0.*chain(:,i)+param_true(i),'--k','LineWidth',2);
grid on; set(gca,'FontSize',20);
title(param_Names{i});
end
%% Density plots
posterior_mode = zeros(N_par,1);
posterior_mean = mean(chain);
figure(2);clf
for i=1:N_par
    subplot(2,2,i); hold on;
    [f,xi] = ksdensity(chain(:,i));
%     [f,xi] = KDE(chain(:,i));
    plot(xi,f,'LineWidth',3);
    plot(0.*xi+param_true(i),linspace(0,max(f),100),'--k','LineWidth',2);
    grid on; set(gca,'FontSize',20);
    title(param_Names{i});
    [~,posterior_id] = max(f);
    posterior_mode(i) = xi(posterior_id);
end

%% Now, evaluate the model of the mode of the posterior distribution
y_mode = call_SEIR([posterior_mode(:)'  infect_cycle],x0,tspace,1:length(tspace));
y_mean = call_SEIR([posterior_mean(:)' infect_cycle],x0,tspace,1:length(tspace));
figure(3);clf;hold on;
plot(tdata,ydata_noisy,'ko','LineWidth',2,'MarkerSize',8);
plot(tspace,ytrue(:,3),'--k','LineWidth',2);
plot(tspace,y_mode,'c','LineWidth',3);
plot(tspace,y_mean,'m','LineWidth',3);
grid on; set(gca,'FontSize',20);
legend('Data','Truth','Mode','Mean');
ylabel('Infected');
xlabel('Time (days)')