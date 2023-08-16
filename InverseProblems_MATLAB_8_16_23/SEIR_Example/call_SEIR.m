% The residual sum of squares for the loggrowth function

function ymodel = call_SEIR(param,x0,tspace,data_ids)

% Pass in the current parameters and solve the ODE system
N = sum(x0);
f = @(t,x) SEIR_model(t,x,param,N);
[~,y] = ode45(f,tspace,x0);

% Now, compare the model predictions to the data
ymodel = y(data_ids,3);
end