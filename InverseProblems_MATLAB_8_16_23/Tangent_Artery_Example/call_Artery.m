% A wrapper file for getting the artery model to provide output at
% specified points
% param - parameters (stiffness, gamma, and P0)
% Dspace - independent variable
% Dref  - reference diameter
% data_ids - the points of the depdendent variable you want to return

function ymodel = call_Artery(param,Dspace,Dref,data_ids)

% Pass in the current parameters and solve the equation
y = Arterial_Model(Dspace,param,Dref);

% Now, grab a subsample of the model response
ymodel = y(data_ids);
end