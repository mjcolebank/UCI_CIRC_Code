% Metropolis hastings algorithm

function [par_chain,s2chain] = MH_Algorithm(inputs)
param0 = inputs.par0;
num_samples = inputs.nsamp;
ssfun = inputs.ss_fun; %should take pars and data
data  = inputs.data;
sig2 = inputs.measurement_noise;
par_low = inputs.par_low;
par_upp = inputs.par_upp;
update_sig = inputs.update_sig;

% For when we updated the variance estimates
a_noiseVar = 0.001; b_noiseVar = 0.001;
n_data = length(data);
% Keep count of rejections
rejout = 0;
acc = 0;
num_par = length(param0);
oldpar = param0;


% default prior function calculates Gaussian sum of squares
% note that uniform prior has no prior function penality
if isfield(inputs,'prior')
    priorfun = @(th,mu,sig) sum(((th-mu)./sig).^2);
    if ~isfield(inputs,'thetamu') || ~isfield(inputs,'thetasig')
        error('Gaussian prior requires mean and standard deviation; exiting.')
    else
        thetamu  = inputs.thetamu;
        thetasig = inputs.thetasig;
    end
else
    priorfun = @(th,mu,sig) 0;
    thetamu = 0;
    thetasig = 0;
end


ss_old    = ssfun(oldpar,data);
prior_old  = priorfun(oldpar,thetamu,thetasig);
par_chain = zeros(num_samples,num_par);
if update_sig
    s2chain = zeros(num_samples,1);
else
    s2chain = [];
end

% We need a proposal distribution: assume Gaussian with covariance
% propotional to the parameter value
% s
qcov = diag((0.05.*param0).^2);
R = chol(qcov);
%%
par_chain(1,:) = param0;
for i=2:num_samples
    q = randn(1,num_par);
    newpar = oldpar + q*R;
    if any(newpar<par_low) || any(newpar>par_upp)
        ss_new = 10^10; % ss
        prior_new = 0;
    else % inside the boundaries
        ss_new = ssfun(newpar,data);
        prior_new = priorfun(newpar,thetamu,thetasig);
    end % inside/outside boundaries

    if ss_new == 10^10 % outside boundaries
        tst = 0;
    else
        tst = exp(-0.5*((ss_new - ss_old)./sig2 + prior_new-prior_old)); % for ss
    end

    if tst <= 0
        accept = 0;
    elseif tst >= 1
        accept = 1; acc = acc + 1;
    elseif tst > rand(1,1)
        accept = 1; acc = acc + 1;
    else
        accept = 0;
    end

    if accept == 1 % accept proposal
        par_chain(i,:) = newpar;
        oldpar = newpar;
        ss_old = ss_new;
        prior_old = prior_new;
    else % reject
        par_chain(i,:) = oldpar;
        rejout = rejout+1;
        %disp('reject')
    end

    % Update error variance
    if update_sig
%                 if i<burnin % sample sigma2 in sampling phase
%                     s2chain(i) = s2chain(i-1);
%                 else
                    s2chain(i) = 1.0/gamrnd(a_noiseVar+0.5*n_data, 1/(b_noiseVar+0.5*ss_old));
%                 end
        sig2 = s2chain(i);
    end

end
disp(rejout./num_samples)
end
%%
