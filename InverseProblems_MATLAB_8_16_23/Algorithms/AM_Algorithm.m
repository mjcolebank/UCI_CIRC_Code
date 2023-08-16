% Adaptive Metropolis

function par_chain = AM_Algorithm(inputs)
param0 = inputs.par0;
num_samples = inputs.nsamp;
ssfun = inputs.ss_fun; %should take pars and data
data  = inputs.data;
sig2 = inputs.measurement_noise;
par_low = inputs.par_low;
par_upp = inputs.par_upp;

if isfield(inputs,'adaptint')
    adapt_int = inputs.adaptint;
else
    adapt_int = 1000;
end

if isfield(inputs,'burnin')
    burnin = inputs.burnin;
else
    burnin = 1000;
end


% Keep count of rejections
rejout = 0;
acc = 0;
last_cov_update = 0;
num_par = length(param0);
oldpar = param0;

% For AM algorithm
covchain = []; meanchain = []; wsum = [];
qcov_adjust = 1e-8; % epsilon adjustment for chain covariance
qcov_scale = 2.4 / sqrt(num_par); % s_d from recursive covariance formula


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

% We need a proposal distribution: assume Gaussian with covariance
% propotional to the parameter values
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

    % Adaptive metropolis
    if i>burnin
        if mod(i, adapt_int) == 0 % we adapt
            [covchain,meanchain,wsum] = covupd(par_chain((last_cov_update+1):i,:),1, ...
                covchain,meanchain,wsum);

            upcov = covchain; % update covariance based on past samples

            [Ra,p] = chol(upcov);
            if p % singular
                % try to blow it
                [Ra,p] = chol(upcov + eye(num_par)*qcov_adjust);
                if p == 0 % choleski decomposition worked
                    % scale R
                    R = Ra * qcov_scale;
                end
            else
                R = Ra * qcov_scale;
            end
            last_cov_update = i;
        end
    end
    
end
disp(rejout./num_samples)
end
%% Covariance update script; copied from Haario and Laine

function [xcov,xmean,wsum,R]=covupd(x,w,oldcov,oldmean,oldwsum,oldR)
%COVUPD covariance update
% [xcov,xmean,wsum]=covupd(x,w,oldcov,oldmean,oldwsum)

% optionally updates also the Cholesky factor R

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:34 $

[n,p]=size(x);
if n == 0 % nothing to update with
    xcov = oldcov; xmean = oldmean; wsum = oldwsum;
    return
end

if nargin<2 | isempty(w)
    w = 1;
end
if length(w) == 1
    w = ones(n,1)*w;
end

if nargin < 6 | isempty(oldR)
    R = [];
else
    R = oldR;
end

if nargin>2 & ~isempty(oldcov) % update

    for i=1:n
        xi     = x(i,:);
        wsum   = w(i);
        xmeann = xi;
        xmean  = oldmean + wsum/(wsum+oldwsum)*(xmeann-oldmean);

        if ~isempty(R)
            R = cholupdate(sqrt((oldwsum-1)/(wsum+oldwsum-1))*R, ...
                (xi-oldmean)'* ...
                sqrt((wsum*oldwsum)/(wsum+oldwsum-1)/(wsum+oldwsum)));
        end

        xcov =  (oldwsum-1)./(wsum+oldwsum-1).*oldcov + ...
            wsum.*oldwsum/(wsum+oldwsum-1)./(wsum+oldwsum) .* ...
            ((xi-oldmean)' *(xi-oldmean));
        wsum    = wsum+oldwsum;
        oldcov  = xcov;
        oldmean = xmean;
        oldwsum = wsum;
    end

else % no update

    wsum  = sum(w);
    xmean = zeros(1,p);
    xcov  = zeros(p,p);
    for i=1:p
        xmean(i) = sum(x(:,i).*w)./wsum;
    end
    if wsum>1
        %%% (wsum-oldwsum/wsum)
        for i=1:p
            for j=1:i
                xcov(i,j) = (x(:,i)-xmean(i))' * ((x(:,j)-xmean(j)).*w)./(wsum-1);
                if (i ~= j)
                    xcov(j,i) = xcov(i,j);
                end
            end
        end
    end

    if nargout>3
        [R,p] = chol(xcov);
        if p~=0
            R=[];
        end
    end

end
end
