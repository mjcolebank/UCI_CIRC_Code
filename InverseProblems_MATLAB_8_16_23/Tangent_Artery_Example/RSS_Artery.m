% The residual sum of squares for the loggrowth function


function J = RSS_Artery(param,data,Dspace,Dref,data_ids)

% Pass in the current parameters and solve the ODE system
ymodel = call_Artery(param,Dspace,Dref,data_ids);

% Calculate the RSS and return the value to the optimizer you are using
if isfield(data,'ydata')
    J = (data.ydata(:)-ymodel(:))'*(data.ydata(:)-ymodel(:));
else
   
J = (data(:)-ymodel(:))'*(data(:)-ymodel(:));
end

end