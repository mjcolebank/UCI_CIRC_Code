% Computed the local sensitivity of the objective function at a specified
% point and return the sensitivities and fisher information matrix

function [sens,F] = local_sensitivity(param,f,ynom,step_size)
% Use a centered difference method for approximating model sensitivity
num_par = length(param);
size_model = length(ynom);
par0       = param;
I_mat      = eye(num_par);
sens       = zeros(size_model,num_par);

% if min(abs(par0))./max(abs(par0))<0.01 % Magnitudes are too different, logscale
%     par_sign = sign(par0);
%     for i=1:num_par
%     step_plus  = log(abs(par0)) + step_size.*I_mat(i,:);
%     step_minus = log(abs(par0)) - step_size.*I_mat(i,:);
%     yplus      = f(exp(step_plus).*par_sign);
%     yminus     = f(exp(step_minus).*par_sign);
%     sens(:,i)  = (yplus-yminus)./(2.*step_size);
%     end
%     sens = sens./par0;
% else
    for i=1:num_par
        step_plus  = par0 + step_size.*I_mat(i,:);
        step_minus = par0 - step_size.*I_mat(i,:);
        yplus      = f(step_plus);
        yminus     = f(step_minus);
        sens(:,i)  = (yplus-yminus)./(2.*step_size);
    end
% end
F = sens'*sens;

end