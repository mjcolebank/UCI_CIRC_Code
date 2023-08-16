% In this example, we will:
% (1) Generate data from a single ODE equation and add WGN
% (2) Construct a "cost function" and perform OLS
% (3) Interpret our parameters after model calibration
clear; clc; close all;
addpath('../Algorithms/')
%% A) First, we will define our equation 
npts  = 101;
Tend = 150;
tspace = linspace(0,Tend,npts); % This will be our time interval in days

% Define our initial conditions  for the ODE model and the parameters
S0 = 999; E0 = 0; I0 = 1; R0 = 0;
beta = 0.8; mu = 0.1; gamma = 0.1; delta = 0.4; infect_cycle = 25; % 

% To make things easier, stack our parameter values into a vector
param = [beta mu gamma delta infect_cycle];
x0 = [S0; E0; I0; R0]; N = sum(x0);

f = @(t,x) SEIR_model(t,x,param,N);
[t,ytrue] = ode45(f,tspace,x0);

%% B) Now, we will generate noisy, subsampled data, from our true model
data_ids = 5:5:npts; tdata = tspace(data_ids); 
ydata_clean = ytrue(data_ids,3);

% Take the clean add Gaussian noise with a fixed error variance
noise_var = 6; noise_mean = 8;
ydata_noisy = ydata_clean+normrnd(noise_mean,noise_var,length(data_ids),1);
% Make sure signal is positive
ydata_noisy = max(ydata_noisy,0);


%% C) Consider model calibration
% NOTE: we don't estimate the infection time here, so just add it as the
% last component of the vector 'q'
RSS_func = @(q) RSS_SEIR([q infect_cycle],ydata_noisy,x0,tspace,data_ids);

% Now, use a built in optmization routine, like fminsearch or fmincon
par_upper = [0.95 0.95 0.95 0.95]; %Upper bound
par_lower = [0.01 0.01 0.01 0.01];  %Lower bound
par_guess = [0.75 0.05 0.15 0.3];   %Initial Guess
[paramopt,RSSopt] = fmincon(RSS_func,par_guess,[],[],[],[],par_lower,par_upper);
%% D) Frequentist confidence intervals
% First, calculate the local sensitivity matrix at the optimal parameter
% vector
N_data = length(data_ids); % Number of data points
N_par  = length(paramopt); % Number of parameters
sens_func = @(q) call_SEIR([q infect_cycle],x0,tspace,data_ids); %Function call for sensitivity analysis
ynom      = sens_func(paramopt);
% Now, compute the sensitivity of the model at the optimal parameter value
% using derivative based local sensitivity analysis
[sens,F] = local_sensitivity(paramopt,sens_func,ynom,1e-3);

% Estimate the unbiased error variance
sigma2_est = RSSopt./(N_data-N_par);

% Finally, get the covariance estimator
cov_est    = sigma2_est.*inv(F);

% Parameter confidence intervals
alpha_conf = 0.975;
par_SE = sqrt(diag(cov_est)).*tinv(alpha_conf,N_data-N_par);
par_CI = [paramopt'-par_SE paramopt'+par_SE];
disp(par_CI)

% We can also get model confidence and prediction intervals
Yopt = call_SEIR([paramopt infect_cycle],x0,tspace,data_ids); % Optimal model prediction
Y_SE_CI = zeros(N_data,1); %Standard error in output space for confidence intervals
Y_SE_PI = zeros(N_data,1); %Standard error in output space for prediction intervals
for i=1:N_data
    Y_SE_CI(i) = sqrt(sens(i,:)*cov_est*sens(i,:)').*tinv(alpha_conf,N_data-N_par);
    Y_SE_PI(i) = sqrt(sigma2_est+sens(i,:)*cov_est*sens(i,:)').*tinv(alpha_conf,N_data-N_par);
end
Y_CI = [Yopt+Y_SE_CI,Yopt-Y_SE_CI]; %Confidence Intervals
Y_PI = [Yopt+Y_SE_PI,Yopt-Y_SE_PI]; %Prediction Intervals (need error variance here)
%%
figure(1); hold on;
h1 = plot(tdata,Yopt,'k','LineWidth',3);
h2 = plot(tdata,Y_CI,'--b','LineWidth',3);
h3 = plot(tdata,Y_PI,'--m','LineWidth',3);
h4 = plot(tdata,ydata_noisy,'ko','LineWidth',3,'MarkerSize',8);
grid on; set(gca,'FontSize',20);
legend([h1 h2(1) h3(1) h4],{'Opt','CI','PI','Data'},'Location','northeast')