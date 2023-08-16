% In this example, we will generate output from an algebraic model
% describing blood vessel mechanics
clear; clc; close all;
%% First, we will define our equation - a vasoactive pressure-area relationship
%       P = stiffness.*tan(pi.*(D/Dref - 1)/gamma);
% based on Qureshi et al. 2019, doi:10.1007/s10237-018-1078-8

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

% Visualize the model output 
P = Arterial_Model(Dspace,param,Dref);

% Pressure-Diameter
figure(1);clf;
plot(Dspace,P,'LineWidth',3);
xlabel('Diameter (cm)')
ylabel('Pressure (mmHg)')
grid on; set(gca,'FontSize',20); axis tight;

% Pressure-Strain
figure(2);clf;
plot(Dstrain,P,'LineWidth',3);
xlabel('Strain (ND)')
ylabel('Pressure (mmHg)')
grid on; set(gca,'FontSize',20); axis tight;

