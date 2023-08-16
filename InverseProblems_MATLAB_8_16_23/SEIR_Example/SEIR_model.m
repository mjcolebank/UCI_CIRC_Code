%% Epidemic Model - differential equation system
% the general equation is an ODE given by compartment modeling
%       dS/dt = S
% where r (1/days) is the growth rate and K (# of organisms) is the carrying
% capacity
% Inputs: x - the current population
%         t - the current time (needed for the ODE solver)
%         param - a vector of parameters (i.e., [r K])
% Outputs: dx - the right hand side of the ODE

function dx = SEIR_model(t,x,param,N)
beta        = param(1);
mu          = param(2);
gamma       = param(3);
delta       = param(4);
infect_cycle = param(5);

if infect_cycle==0
    delta_f = delta;
else
    delta_f = max(delta.*sin(t./infect_cycle),0);
end

S = x(1);
E = x(2);
I = x(3);
R = x(4);

dx = zeros(4,length(t));
dx(1) = -beta.*S.*I./N - mu.*S + mu.*N;
dx(2) =  beta.*S.*I./N - (mu+delta_f).*E;
dx(3) =  delta_f.*E - (gamma+mu).*I;
dx(4) =  gamma.*I - mu.*R;
end