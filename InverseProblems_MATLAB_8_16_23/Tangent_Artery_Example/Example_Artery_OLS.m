% In this example, we will:
% (1) Generate data from a single equation and add WGN
% (2) Construct a "cost function" and perform OLS
% (3) Interpret our parameters after model calibration
clear; clc; close all;
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

% Take the data and add Gaussian noise with a fixed error variance
noise_var = 2;
Pdata_noisy = Pdata_clean+normrnd(0,noise_var,1,length(data_ids));

% Visualize the noisy data on top of the model
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
disp([param; paramopt]') % Display estimated parameters vs true parameters

% Look at the model predictions at the optimal parameter vector
Popt = Arterial_Model(Dspace,paramopt,Dref);

figure(1); hold on;
plot(Dspace,Popt,'c','LineWidth',3);
grid on; set(gca,'FontSize',20);
ylabel('Pressure (mmHg)'); xlabel('Diameter (cm)')

% Always check your residuals!
figure(2); hold on;
plot(Ddata,Popt(data_ids)-Pdata_noisy,'ko','LineWidth',3,'MarkerSize',12);
grid on; set(gca,'FontSize',20);
ylabel('Residual'); xlabel('Diameter (cm)')