% In this example, we will:
% (1) Generate data from a single ODE equation and add WGN
% (2) Construct a "cost function" and perform OLS
% (3) Interpret our parameters after model calibration
clear; clc; close all;
addpath('../Algorithms/')
%% A) First, we will define our equation - the Logistic Growth
npts   = 101; tspace = linspace(0,100,npts); % This will be our time interval in days

% Define our initial conditions  for the ODE model and the parameters
x0 = 2;  r = 0.12; K = 150;
param = [r K];
f = @(t,x) logistic_growth_model(t,x,param);
[t,y] = ode45(f,tspace,x0);

%% B) Now, we will generate noisy, subsampled data, from our true model
data_ids = 5:5:npts; tdata = tspace(data_ids); ydata_clean = y(data_ids);

% Take the clean add Gaussian noise with a fixed error variance
noise_var = 5; ydata_noisy = ydata_clean+normrnd(0,noise_var,length(data_ids),1);

%% C) Consider model calibration
% Now, use a built in optmization routine, like fminsearch or fmincon
RSS_func = @(q) RSS_loggrowth(q,ydata_noisy,x0,tspace,data_ids);
par_upper = [0.8 300]; %Upper bound
par_lower = [0.01 10];  %Lower bound
par_guess = [0.1 100];   %Initial Guess
[paramopt,RSSopt] = fmincon(RSS_func,par_guess,[],[],[],[],par_lower,par_upper);
disp([param; paramopt]') % Display estimated parameters vs true parameters
%% D) Frequentist confidence intervals
% First, calculate the local sensitivity matrix at the optimal parameter
% vector\
N_data = length(data_ids); % Number of data points
N_par  = length(param);    % Number of parameters
sens_func = @(q) call_loggrowth(q,tspace,x0,data_ids); %Function call for sensitivity analysis
ynom      = sens_func(paramopt); % Solution at the optimal parameter vector

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
Yopt = call_loggrowth(paramopt,tspace,x0,data_ids); % Optimal model solution
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
legend([h1 h2(1) h3(1) h4],{'Opt','CI','PI','Data'},'Location','southeast')