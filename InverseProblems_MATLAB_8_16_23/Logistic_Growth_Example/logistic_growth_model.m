%% Logistic growth model - differential equation system
% the general equation is an ODE given by
%       dx/dt = r*x*(1-x/K)
% where r (1/days) is the growth rate and K (# of organisms) is the carrying
% capacity
% Inputs: x - the current population
%         t - the current time (needed for the ODE solver)
%         param - a vector of parameters (i.e., [r K])
% Outputs: dx - the right hand side of the ODE

function dx = logistic_growth_model(t,x,param)

% Version 1: growth with two parameters
r = param(1);
K = param(2);

dx = r.*x.*(1-x./K);

end
% Version 2: Growth with oscillatory harvesting cycles
% r = param(1);
% K = param(2);
% M = param(3);
% 
% dx = r.*x.*(1-x./K) - x.*M.*sin(t./20);
% end