% Call the logistic growth model and return the ODE solution

function ymodel = call_loggrowth(param,tspace,x0,data_ids)
% Pass in the current parameters and solve the ODE system
f = @(t,x) logistic_growth_model(t,x,param);
[~,y] = ode45(f,tspace,x0);

% Now, get model predictions where there is data
ymodel = y(data_ids);
end