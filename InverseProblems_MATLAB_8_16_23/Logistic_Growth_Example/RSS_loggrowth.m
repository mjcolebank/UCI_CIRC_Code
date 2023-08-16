% The residual sum of squares for the loggrowth function


function J = RSS_loggrowth(param,data,x0,tspace,data_ids)

% Pass in the current parameters and solve the ODE system
ymodel = call_loggrowth(param,tspace,x0,data_ids);

% Calculate the RSS and return the value to the optimizer you are using
J = (data(:)-ymodel(:))'*(data(:)-ymodel(:));

end