% In this example, we will:
% (1) Generate data from a single ODE equation and add WGN
% (2) Construct a "cost function" and perform OLS
% (3) Interpret our parameters after model calibration
% (4) Quantify uncertainty using frequentist statistical theory
clear; clc; close all;
addpath('../Algorithms/')
%% A) First, we will define our equation - the Logistic Growth
% Define the abscissa for the model (here, Diameter of a blood vessel)
npts   = 101;
Dref   = 6.00;
Dmax   = 7.00;

Dspace = linspace(Dref,Dmax,npts); % This will be x values
Dstrain = (Dspace-Dref)./Dref;     % This is the strain, useful for biomechanics folks

% Define our parameters
stiff = 30;   % Stiffness of the artery (mmHg)
gamma = 0.5;  % Inflection point in the tangent curve (dimensionless)
P0    = 80;   % Reference pressure at D=Dref (mmHg)

% To make things easier, stack our parameter values into a vector
param = [stiff gamma P0];
P = Arterial_Model(Dspace,param,Dref);
%% B) Now, we will generate noisy, subsampled data, from our true model
data_ids = 1:10:npts;
Ddata = Dspace(data_ids);
Pdata_clean = P(data_ids);

% Take the clean add Gaussian noise with a fixed error variance
noise_var = 2;
Pdata_noisy = Pdata_clean+normrnd(0,noise_var,1,length(data_ids));

% Visual the noisy data on top of the model
figure(1);clf;hold on;
plot(Dspace,P,'k','LineWidth',3);
plot(Ddata,Pdata_clean,'ko','LineWidth',3,'MarkerSize',12); 
plot(Ddata,Pdata_noisy,'--r+','LineWidth',3,'MarkerSize',12); 
ylabel('Pressure (mmHg)'); xlabel('Diameter (cm)');
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
RSS_func = @(q) RSS_Artery(q,Pdata_noisy,Dspace,Dref,data_ids);

% Now, use a built in optmization routine, like fminsearch or fmincon
par_upper = [70 1  120]; %Upper bound
par_lower = [20 0.01 50];  %Lower bound
par_guess = [20 0.9  60];   %Initial Guess
[paramopt,RSSopt] = fmincon(RSS_func,par_guess,[],[],[],[],par_lower,par_upper);
%% D) Frequentist confidence intervals
% First, calculate the local sensitivity matrix at the optimal parameter
% vector
N_data = length(data_ids); % Number of data points
N_par  = length(param);    % Number of parameters
sens_func = @(q) call_Artery(q,Dspace,Dref,data_ids); %Function call for sensitivity analysis
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
Yopt = call_Artery(paramopt,Dspace,Dref,data_ids)'; % Optimal model prediction
Y_SE_CI = zeros(N_data,1); %Standard error in output space for confidence intervals
Y_SE_PI = zeros(N_data,1); %Standard error in output space for prediction intervals
for i=1:N_data
    Y_SE_CI(i) = sqrt(sens(i,:)*cov_est*sens(i,:)').*tinv(alpha_conf,N_data-N_par);
    Y_SE_PI(i) = sqrt(sigma2_est+sens(i,:)*cov_est*sens(i,:)').*tinv(alpha_conf,N_data-N_par);
end
Y_CI = [Yopt+Y_SE_CI,Yopt-Y_SE_CI]; %Confidence Intervals
Y_PI = [Yopt+Y_SE_PI,Yopt-Y_SE_PI]; %Prediction Intervals (need error variance here)
%%
figure(2); hold on;
h1 = plot(Ddata,Yopt,'k','LineWidth',3);
h2 = plot(Ddata,Y_CI,'--b','LineWidth',3);
h3 = plot(Ddata,Y_PI,'--m','LineWidth',3);
h4 = plot(Ddata,Pdata_noisy,'ko','LineWidth',3,'MarkerSize',8);
grid on; set(gca,'FontSize',20);
legend([h1 h2(1) h3(1) h4],{'Opt','CI','PI','Data'},'Location','southeast')