function [x,P,resids,L,Scov] = C19filt_kalman(KF,y,x0,P0,beta0)
%C19FILT_KALMAN Kalman filter with state dependent covariance.
%   [x,P,resids] = C19FILT_KALMAN(KF,y,x0,P0) uses the Kalman filter
%   defined by KF and the measurements defined in y, to produce state
%   estimates x and covariance matrices P. x0 and P0 set the initial
%   state and covariance, respectively. The output resids is the
%   residuals between data and prior.
%
%   [x,P,resids,L] = C19FILT_KALMAN(KF,y,x0,P0) additionally computes the
%   marginal likelihood L of data y given the model (x0,P0). L is a
%   length size(y,2) vector with L(m) the marginal likelihood
%   contribution of y(:,m).
%
%   [x,P,resids,L,Scov] = C19FILT_KALMAN(KF,y,x0,P0) also returns the
%   innovation covariance Scov for each time index, i.e. Scov(:,:,m) is
%   the covariance for index m.
%
%   [x,P] = C19FILT_KALMAN(KF,y,x0,P0,beta0) same as case 1, but with
%   dynamical beta.
%
%   The state estimate x has the following structure: x(:,m,1) and
%   x(:,m,2) contain the a priori- and a posteriori state estimates at
%   index m, from the measured data y(:,1:m-1) and y(:,1:m),
%   respectively. P(:,:,m,1) and P(:,:,m,2) contain the a priori- and
%   a posteriori covariances correspondingly.
%
%   The Kalman filter is defined by the following fields:
%
%   Field   Definition
%   ------------------------------------------------------------------
%   F       State update matrix, x_{i+1} = F*x_i.
%
%   H       Measurement operator, y_i = H*x_i.
%
%   R0      The measurement noise covariance computed by the
%   rdiag   formula R = R0+rdiag*(H*xx).^2.
%
%   Qp      Local contribution of the covariance from (local) state
%           x_i according to Q_i = Qp(:,:,i)*x_i, Q_local =
%           sum(Q_i).
%
%   Q0      Q0, qdiag, and Qv define additional uncertainties
%   qdiag   according to Q = Q + Q0 + qdiag*diag(Qv*xx.^2).
%   Qv
%
%   See also C19FILT.

% Hidden support for exception log:
%
%   clear C19filt_kalman; % clears the log
%   [...] = C19filt_kalman(...); % use filter
%   exclog = C19filt_kalman; % retrieves the log
%
% The log counts the different types of exceptions and also records
% the index in time when the exception was recorded.

% S. Engblom 2021-05-15 (exception log)
% R. Eriksson 2021-01-28 (merge with dynamical_beta branch)
% H. Runvik 2020-11-19 (New qdiag formula)
% S. Engblom 2020-10-17 (Bugfix prior on (x0,P0))
% S. Engblom 2020-10-16 (Additional output L)
% S. Engblom 2020-10-15 (Speed improvements)
% H. Runvik 2020-09-14

if nargin < 5 || isempty(beta0)
  beta = [];
  betaIdx = [];
else
  beta = beta0.beta;
  betaIdx = beta0.betaIdx;
end

try 
 exceptioncheck = isfield(KF,'exception') && nargout > 3;
catch % allows for hidden support.
 exceptioncheck = false;
end

earlylog = nargin == 0 && nargout == 1; % return the exceptions log
if exceptioncheck || earlylog
 % exception log:
 persistent exclog;
 if earlylog
  % special syntax to retreieve the log:
  x = exclog;
  return;
 end
 
 exc = KF.exception; % convenience
 SAFETY = double(realmin('single'))*realmax;
 % (like "single precision away from being double precision -inf")
end

% Kalman filter: dynamics and measurements
F = KF.F;
% (change into more effective format:)
Qp = sparse(reshape(KF.Qp,[],size(KF.Qp,3)));
Q0 = KF.Q0;

% sizes
nstate = size(KF.Qp,3);
nnode = size(F,1)/nstate;
ntime = size(y,2);
Id = speye(size(F));

% *** more effective format possible (sparse)?
Qv = kron(eye(nnode),KF.Qv);

% build block indexes (iiQ,jjQ)
Q_ = cell(1,nnode);
Q_(:) = {sparse(ones(nstate))};
Q_ = blkdiag(Q_{:});
[iiQ,jjQ,~] = find(Q_);
% *** [iiQ,jjQ] should be built elsewhere?

H = KF.H;
R0 = KF.R0;
Rdiag = sparse(1:size(H,1),1:size(H,1),KF.rdiag);

% allocate filter output
x = zeros(size(F,1),ntime,2);
P = zeros(size(F,1),size(F,2),ntime,2);
resids = zeros(size(y,1),ntime);
Scov=zeros(size(H,1),size(H,1),ntime);

% initial state
xx = x0;
PP = P0;
L = [];
if nargout > 3
  L = zeros(1,size(y,2));
  d2pi = size(y,1)*log(2*pi);
end

for m = 1:ntime
  if ~isempty(beta)
    F(betaIdx) = beta(m);
  end
  % prior on states *after* the initial state
  if m > 1
    % additive noise model (covariance)
    Q = Qp*reshape(xx,nstate,[]); % effective covariance operator
    Q = sparse(iiQ,jjQ,Q,size(F,1),size(F,2));

    % add offset and rate-dependent covariance
    Q = Q+Q0+KF.qdiag*diag(Qv*xx.^2);

    % a priori
    xx = F*xx;
    PP = F*PP*F'+Q; % [profile: shared 2nd worst]
  end

  % a priori
  x(:,m,1) = xx;
  P(:,:,m,1) = PP;

  % exceptional: check for NaN in measurement
  yNum = y(:,m);
  Hmod = H;
  if any(isnan(yNum))
    ixn = find(~isnan(yNum));
    Z = sparse(ixn,ixn,1,numel(yNum),numel(yNum));
    % remove them:
    yNum = Z*yNum;
    Hmod = Z*H;
  end

  % pre-fit residual
  yy = yNum-Hmod*xx;
  resids(:,m) = yy;

  % update
  R = R0+Rdiag.*(Hmod*xx).^2;
  S = Hmod*PP*Hmod'+R;
  K = PP*Hmod'/S; % [profile: 3rd worst]
  
  Scov(:,:,m) = S;
  
  % updated likelihood
  if ~isempty(L)
    % (note: unclear what to return if det < 0, but likely better to have
    % it non-silent and investigate...)
    L(m) = -0.5*(yy'*(S\yy)+log(det(S))+d2pi);
  end

  % a posteriori
  xx = xx+K*yy;
  PP = (Id-K*Hmod)*PP; % [profile: 1st worst]
  x(:,m,2) = xx;
  P(:,:,m,2) = PP;
  
  if exceptioncheck
    % bad L or bad state (xx,PP)?
    if L(m) < -SAFETY || ... % very small likelihood?
          ~isreal(L(m)) || ...       % oops! went imaginary!
          isnan(L(m)) || ...         % oops! went NaN!
          any(xx < exc.LB) || ...    % mean outside bounds [LB,UB]
          any(exc.UB < xx) || ... 
          any(diag(PP) < 0) || ...   % no negative diagonals
          any(exc.AbsMagn < xx & xx.^2 < exc.SDFAC^2*diag(PP))
      % final check is: *above* AbsMag *and* large portion of probability
      % mass is negative (norminv(0.4,0,1) = -0.2533 so ~40%
      % probability is negative with SDFAC = 0.25)
      if isempty(exclog)
        % allocate
        exclog.count = zeros(1,7);
        exclog.index = zeros(1,0);
        exclog.state = zeros(numel(xx),0);
      end
      % log the exception
      if L(m) < -SAFETY, exclog.count(1) = exclog.count(1)+1; end
      if ~isreal(L(m)), exclog.count(2) = exclog.count(2)+1; end
      if isnan(L(m)), exclog.count(3) = exclog.count(3)+1; end
      if any(xx < exc.LB), exclog.count(4) = exclog.count(4)+1; end
      if any(exc.UB < xx), exclog.count(5) = exclog.count(5)+1; end
      if any(diag(PP) < 0), exclog.count(6) = exclog.count(6)+1; end
      if any(exc.AbsMagn < xx & xx.^2 < exc.SDFAC^2*diag(PP))
        exclog.count(7) = exclog.count(7)+1;
      end
      exclog.index = [exclog.index m];
      exclog.state = [exclog.state xx];
      exclog.exc = exc;
      % then throw!
      L(m:end) = -inf;
      return;
    end
  end
end
