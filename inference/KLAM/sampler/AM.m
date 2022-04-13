function [thetas, sl, param, outverb] = AM(nMCMC, thetas, sl, param, funcs, Ydata, start, stop)
%AM finds sampels from the (approximate) posterior using Adaptive
%Metropolis.
%   [THETAS, SL] = AM(NMCMC, THETAS, SL, PARAM, LOGKL, PRIOR, START, STOP)
%   returns samples from the approximate posterior using the KLAM method.
%
%   --- output ---
%   THETAS - samples from the posterior. Remember to check for burnin
%   period.
%   KL - (Synthetic) likelihood for each sample in THETAS
%   OUTVERB - diagnostic variables that are also displayed in prompter if
%   PARAM.verb = true.
%
%   --- input ---
%   NMCMC = the length of the chain.
%   THETAS = matrix with current chain, can be empty except for first
%          element.
%   KL = synthetic likelihood for the current chain.
%   PARAM = AM method parameters.
%   funcs = {LOGKL -function to compute the synthetic log likelihood.
%            PRIOR = the prior density.}
%   YDATA = data for inference.
%   START = where to start in the given chain.
%   STOP  = where to stop the the given chain.
%
%

% R. Eriksson 201006
if nargin < 8 % starting up a new chain, otherwise, these have to be set.
  start = 2;
  stop = nMCMC;
end



% starting guess for AM
sigma = param.sigma;
xcov = param.xcov;
xbar = param.xbar;

% initial sample
theta = thetas(:,start-1)';

% prior density of initial sample
prior_ = param.prior;
prior_old = prior_.sumlog;

% likelihood of initial sample
L = sl(:,start-1);

% debug information.
accept_rate    = 0;
progress       = 0;
priorfail      = 0;
earlyReject    = 0;
naccept        = 0;
nreject        = 0;

% keeping track of early rejects
part1 = round((stop-start)/10); % first 10%
part2 = round((stop-start)*9/10); % 10-90%
part3 = stop-start; % final 10%
for i = start:stop

  %% AM - proposal step
  % Are we post the burn-in period? Then update the covariance matrix
  % Else, keep the inital proposal covariance.
  if i > param.i0
    % recursively update the covariance and mean.
    [xcov, xbar] = covRec(thetas(:,i-1)',xcov,xbar,i-1);
    sigma = param.Sd*xcov + param.Sd*param.e*eye(numel(theta));
  end

  % stepping, new samples are generated by adding a MVN random variable
  % which is centered at 0 with covariance sigma.
  eps = mvnrnd(zeros(length(theta),1), sigma);
  theta_new = theta + eps;

  % c19filt used a struct format, therefore, construct the struct using
  % mat2stuct with the new proposal (in matrix form).

  theta_new_struct = funcs.m2s(theta_new');

  % check the prior density of the proposal. If it's 0 or log(prior) =
  % -Inf, then reject before computing the likelihood, as it's outside of
  % the support of the prior.
  prior = funcs.prior(theta_new_struct);
  prior_new = getfield(prior,'sumlog');



  % check if prior is inside of support
  if ~isinf(prior_new)
    % Find the likelihood of the current slab.
    [~,~,L_new] = funcs.logKL(theta_new_struct, Ydata);
    % compute the cumulative likelihood "so far".
    L_sum = sum(L_new);

    % For some rate combinations, the likelihood becomes invalid, i.e.,
    % imaginary or nan.
    if abs(imag(L_sum)) == 0 && ~isnan(L_sum) && ~isinf(L_sum)
      % Compute the acceptance probability. which is
      % min{1, exp( [L_new + prior_new] - [L_old + prior_old] ) }

      %
      u = rand;
      P_new = L_sum + prior_new;
      P_old = sum(L) + prior_old;

      % finally compute the acceptance probability
      a = min(1, exp(P_new - P_old));

      keep = u < a; % accept
    else
      earlyReject = earlyReject + 1;
      keep = false;
    end % if abs(imag(...


    if keep
      theta = theta_new;
      L = L_new;
      prior_old = prior_new;
      prior_ = prior;
      naccept = naccept+1;
    end % if keep
  else % is inf
    priorfail = priorfail + 1;
    keep = false;
  end % end exp


  if ~keep
    nreject = nreject+1;
  end

  % store the markov chains.
  sl(:,i) = L;
  thetas(:,i) = theta;




  % For futher debugging, we can display the current acceptance rate level
  % and the progress.
  if mod(naccept+nreject,1e2) == 0
    accept_rate = naccept/(naccept+nreject);
    if param.verb
      progress = (i-start)/(stop-start);
      disp([accept_rate, progress, param.Sd])
    end
  end

  % back-up of the current chain if something breaks.
  if mod(i,1e4) == 0
    save([param.savepath 'AM_temp.mat'])
  end

  % store for later analyzis
  if i == part1
    earlyReject1 = earlyReject;
    nreject1 = nreject;
  elseif i == part2
    earlyReject2 = earlyReject;
    nreject2 = nreject;
  elseif i == part3
    earlyReject3 = earlyReject;
    nreject3 = nreject;
  end
end


% to continue the chain, say we run for 1e5 values, and then want 1e5
% values more. We keep some method parameters that the final sample had.
param.xcov = xcov;
param.xbar = xbar;
param.prior = prior_;

% filling debugging values.
outverb = struct();
outverb.accept_rate = naccept/(naccept+nreject);
outverb.prior_reject_rate = priorfail/(naccept+nreject);

try
  % early rejection rate based on the three sections.
  if part1 > 0 && part2 > 0 && part3 > 0
    outverb.early_reject_rate_10 = earlyReject1 / nreject1;
    outverb.early_reject_rate_10_90 = ...
      (earlyReject2 - nreject1*outverb.early_reject_rate_10)...
      /(nreject2-nreject1);
    outverb.early_reject_rate_90_100 = ...
      (earlyReject3 - ...
      (nreject2-nreject1)*outverb.early_reject_rate_10_90 - ...
      nreject1*outverb.early_reject_rate_10)/(nreject3-nreject2);
  end

  outverb.early_reject_rate = earlyReject/nreject;
catch
  outverb.early_reject_rate_10 = [];
  outverb.early_reject_rate_10_90 = [];
  outverb.early_reject_rate_90_100 = [];
  outverb.early_reject_rate = [];
end
if param.verb
  disp('done')
  disp(['accept_rate ' num2str(outverb.accept_rate)]);
  disp(['priorfailed ' num2str(outverb.prior_reject_rate)]);
  disp(['earlyReject ' num2str(outverb.early_reject_rate)]);
end

end