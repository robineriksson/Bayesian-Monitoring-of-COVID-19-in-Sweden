function [xx,PP] = initC19filt(F,H,y,CStates,Y0,Y0Cov,PScales)
%INITC19FILT Initial state estimation for Kalman Covid-19 filter.
%   xx = INITC19FILT(F,H,y,CStates,[],[],PScales) estimates an initial
%   state xx and a corresponding covariance matrix PP for data y, state
%   update matrix F, observation matrix H and covariance parameters
%   PScales. The logical array CStates defines which states in F are
%   cumulative. The vector of measurements y relates to the states as y =
%   H*x.
%
%   The first step of the algorithm removes all cumulative states from
%   F, H and y to produce F_red, H_red and y_red. The dominating
%   eigenvector of F_red is then projected onto the subspace defined
%   by H_red*x_red = y_red, under a positivity constraint. For
%   cumulative states that are measured, the initial state is taken
%   directly as the measurement, and remaining cumulative states are
%   set to zero. By initializing the filter like this, transients are
%   reduced. In particular, if the system would be initialized from
%   the eigenvector and simulated without perturbations, the relative
%   magnitude of the non-cumulative states would remain constant.
%
%   The covariance matrix PP is assumed to be diagonal with entries
%   proportional to the states squared. The optional input length 2
%   vector PScales contains the values [observed unobserved] (default
%   [1 10]), which defines the scaling of the covariances
%   corresponding to states that are directly measured and,
%   respectively, not measured.
%
%   The user can also supply the initialization with Y0 and Y0Cov whichs
%   holds information from a previously run filter, and the user wants to
%   pick up where they left off.
%
%   Assumptions:
%
%   -All modes of the system dynamics other than the dominating one is
%   neglected (which would be the case if the system has been running
%   without perturbations for a sufficently long period).
%
%   -The measurements are accurate, i.e., the elements of xx will
%   satisfy H*xx = y.
%
%   Known issue:
%
%   -If the dominating eigenvalue is less than one, some
%   parametrizations of the F matrix results in the dominating
%   eigenvector being confined to the H-W-P-subspace (occurs when one
%   of these states is the "slowest to empty"). This will give
%   unreasonable zero initialization of I, E, A and phi.

% H. Runvik 2021-03-04 (Covariance init)
% H. Runvik 2021-02-22 (NaN handling)
% H. Runvik 2020-10-25 (Revision)
% H. Runvik 2020-10-13 (Revision)
% H. Runvik 2020-10-02 (Revision)
% H. Runvik 2020-09-14

switch nargin
  case 4
    Y0 = [];
    Y0Cov = [];
    PScales = [0.1 1];
  case 5
    Y0Cov = [];
    PScales = [0.1 1];
  case 6
    PScales = [0.1 1];
end

if ~isempty(Y0) && ~isempty(Y0Cov) % start from provided initialization
  %warning('not yet fully defined');
  
  % safety check, are they reasonble?
  ok = 1;
  % are they of correct size?
  ok = ok && isequal(size(Y0Cov),size(F));
  ok = ok && isequal(numel(Y0), size(F,1));
  
  % is the starting guess in the "ball park" if the data?
  ok = ok && norm(H*Y0 - y) < 1e-2*norm(y);
  if ok
    xx = Y0;
    PP = Y0Cov;
  else
    error('Look at initalization');
  end
else % start from scracth
  % remove cumulative states from F
  F_red = F(~CStates,~CStates);
  
  % remove measurements of cumulative states (assuming directly measured,
  % i.e. not combined with other states)
  H_red = H(:,~CStates);
  H_idx = all(H_red == 0,2);
  DStates = mod(find(H(H_idx,:)')-1,size(H,2))+1;
  
  % remove NaN measurement
  NaN_idx = isnan(y);
  y_red = y(~H_idx&~NaN_idx);
  H_red = H_red(~H_idx&~NaN_idx,:);
  
  % find dominating eigenvalue of F
  [V,D] = eig(F_red);
  aDiag = real(diag(D)); % ignore imaginary part for now
  [~,imax] = max(aDiag);
  Vmax = V(:,imax);
  
  % rescale the dominating eigenvector to fit measured data
  c = H_red*Vmax\y_red;
  Vscaled = Vmax*c;
  
  % Estimated state satisfying H*xx = y, this determines xx up to a
  % vector in the nullspace of H. The remaining vector is determined so
  % that the projections of Vscaled and xx onto the nullspace coincide.
  
  % Unconstrained version:
  H_null = null(full(H_red))';
  S = [H_red; H_null];
  xx_red = S\[y_red; H_null*Vscaled];
  
  % With positivity constraint:
  %opts = optimoptions(@lsqlin,'Display','off');
  %xx_red = lsqlin(eye(size(F_red)),Vscaled,[],[],H_red,y_red,zeros(size(F_red,2),1),[],[],opts);
  
  % Add remaining elements of xx
  xx = zeros(numel(CStates),1);
  xx(~CStates) = xx_red;
  xx(DStates) = y(H_idx);
  
  % check violation of non-complex assumption
  ix = imag(xx);
  rx = real(xx);
  if any(abs(ix) > 0.01*abs(rx))
    % (note: in case this happens we need to think again)
    warning('Large imaginary part detected.');
  end
  xx = rx;
  
  % covariance init
  H_p = H(~NaN_idx,:);
  observed_idx = sum(H_p,1) >= 1;
  pp = zeros(1,size(H,2));
  pp(observed_idx) = PScales(1);
  pp(~observed_idx) = PScales(2);
  PP = diag(xx).^2.*pp;
end
