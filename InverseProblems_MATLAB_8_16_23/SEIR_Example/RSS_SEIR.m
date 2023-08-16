% The residual sum of squares for the loggrowth function

function J = RSS_SEIR(param,data,x0,tspace,data_ids)

% Pass in the current parameters and solve the ODE system
% compare the model predictions to the data
ymodel = call_SEIR(param,x0,tspace,data_ids);

% Calculate the RSS and return the value to the optimizer you are using
J = (data(:)-ymodel(:))'*(data(:)-ymodel(:));

end