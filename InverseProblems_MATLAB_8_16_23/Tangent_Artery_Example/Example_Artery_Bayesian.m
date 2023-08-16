% Try running Metropolis Hastings with the Arterial wall model
clear; clc; %close all;
addpath('../Algorithms/')
param_Names = {'Stiffness','\gamma','P0'};
%% A) First, we will define our equation
% Define the abscissa for the model (here, Diameter of a blood vessel)
npts   = 101; Dref   = 6.00; Dmax   = 7.00;
Dspace = linspace(Dref,Dmax,npts); % This will be x values
% Define our parameters
stiff = 30;   % Stiffness of the artery (mmHg)
gamma = 0.5;  % Inflection point in the tangent curve (dimensionless)
P0    = 80;   % Reference pressure at D=Dref (mmHg)

param_true = [stiff gamma P0];
Ptrue = Arterial_Model(Dspace,param_true,Dref);
%% B) Now, we will generate noisy, subsampled data, from our true model
data_ids = 1:5:npts; Ddata = Dspace(data_ids); Pdata_clean = Ptrue(data_ids);

% Take the clean add Gaussian noise with a fixed error variance
noise_var = 1; Pdata_noisy = Pdata_clean+normrnd(0,noise_var,1,length(data_ids));

%% C) Consider model calibration
% Create an in-line function handle
par_upper = [60 1  120]; %Upper bound
par_lower = [20 0.01 20];  %Lower bound
par_guess = [20 0.45  90];   %Initial Guess
RSS_func = @(q) RSS_Artery(q,Pdata_noisy,Dspace,Dref,data_ids);

% If feasible, optimize once first
% [par_guess,RSSopt] = fmincon(RSS_func,par_guess,[],[],[],[],par_lower,par_upper);

N_par = length(par_guess);
N_data = length(Pdata_noisy);
%%
burnin = 1000;
noise_est = RSS_Artery(par_guess,Pdata_noisy,Dspace,Dref,data_ids)./(N_data-N_par);
% noise_est = 5000;
inputs.par0 = par_guess;
inputs.nsamp = 10000;
inputs.ss_fun = @(q,data) RSS_Artery(q,data,Dspace,Dref,data_ids);
inputs.data   = Pdata_noisy;
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
P_mode = Arterial_Model(Dspace,posterior_mode,Dref);
P_mean = Arterial_Model(Dspace,posterior_mean,Dref);
figure(3);clf;hold on;
plot(Ddata,Pdata_noisy,'ko','LineWidth',2,'MarkerSize',8);
plot(Dspace,Ptrue,'--k','LineWidth',2);
plot(Dspace,P_mode,'c','LineWidth',3);
plot(Dspace,P_mean,'m','LineWidth',3);
grid on; set(gca,'FontSize',20);
legend('Data','Truth','Mode','Mean');
ylabel('Pressure (mmHg)');
xlabel('Diameter (cm)')








