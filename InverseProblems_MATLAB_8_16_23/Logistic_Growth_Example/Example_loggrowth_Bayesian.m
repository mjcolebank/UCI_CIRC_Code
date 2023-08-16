% Try running Metropolis Hastings with the logistic growth model
clear; clc; close all;
addpath('../Algorithms/')
param_Names = {'r','K'};
%% A) First, we will define our equation
% Define the abscissa for the model
npts = 101;
tspace = linspace(0,100,npts);
x0 = 2;
param_true = [0.12 150];
[t,ytrue] = ode45(@logistic_growth_model,tspace,x0,[],param_true);

%% B) Now, we will generate noisy, subsampled data, from our true model
data_ids = 5:5:npts;
tdata    = tspace(data_ids);
noise_var = 2;
ydata_noisy = ytrue(data_ids)+normrnd(0,noise_var,length(data_ids),1);

RSS_func = @(q) RSS_loggrowth(q,ydata_noisy,x0,tspace,data_ids);
par_upper = [0.8 300]; %Upper bound
par_lower = [0.01 10];  %Lower bound
par_guess = [0.1 100];   %Initial Guess

% If feasible, optimize once first
% [paramopt,RSSopt] = fmincon(RSS_func,par_guess,[],[],[],[],par_lower,par_upper);
% [paramopt,RSSopt,~,~,~,~,Hes] = fmincon(RSS_func,par_guess,[],[],[],[],par_lower,par_upper);


N_par = length(par_guess);
N_data = length(ydata_noisy);
%%
burnin = 1000;
noise_est = RSS_loggrowth(par_guess,ydata_noisy,x0,tspace,data_ids)./(N_data-N_par);
% noise_est = 5000; 
inputs.par0 = par_guess;
inputs.nsamp = 10000;
inputs.ss_fun = @(q,data) RSS_loggrowth(q,data,x0,tspace,data_ids);
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

%%
figure(1);clf;
for i=1:N_par
subplot(1,2,i); hold on;
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
    subplot(1,2,i); hold on;
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
y_mode = call_loggrowth(posterior_mode,tspace,x0,1:length(tspace));
y_mean = call_loggrowth(posterior_mean,tspace,x0,1:length(tspace));
figure(3);clf;hold on;
plot(tdata,ydata_noisy,'ko','LineWidth',2,'MarkerSize',8);
plot(tspace,ytrue,'--k','LineWidth',2);
plot(tspace,y_mode,'c','LineWidth',3);
plot(tspace,y_mean,'m','LineWidth',3);
grid on; set(gca,'FontSize',20);
legend('Data','Truth','Mode','Mean');
ylabel('Population');
xlabel('Time (days)')





