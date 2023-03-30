function [Z,covZ,resids,L,S] = C19filt(Y0,Y0Cov,G,rates,D,obsrates,Ydata, ...
    slab,ntime,beta,Q,exception)
%C19FILT Kalman filter for the EngEr Covid-19 model given data.
%   [Z,covZ] = C19FILT(Y0,Y0Cov,G,RATES,D,OBSRATES,YDATA,SLAB,NTIME, ...
%     BETA,Q,EXCEPTION)
%   estimates output states Z = G*Y (see below) given an initial state
%   (Y0,CovY0), model rates RATES, spatial transmission matrix D,
%   observation operator given by OBSRATES, and for data YDATA over a
%   window of time 1:NTIME. Including EXCEPTION allows for tuning of the
%   exceptions used by the Kalman filter for early rejection. If left
%   empty, standard exceptions are deployed.
%
%   If (Y0,CovY0) are both empty, INITC19FILT is used to estimate an
%   initial state. RATES and OBSRATES are descriped in GETC19SYST and
%   GETC19OBS, respectively. D is an NNODES-by-NNODES sparse matrix
%   describing the transmission of infectious pressure phi between
%   nodes. The last four inputs are optional. SLAB defines the
%   division of the time window into slabs; its elements define which
%   slab the corresponding point in time belongs to. If SLAB is not
%   set, a single slab is used. The NTIME parameter allows the time
%   window to exceed the available data and can be used for prediction
%   experiments. BETA is a struct holding fields {beta, betaIdx} for
%   supporting cases of dynamic/daily beta. Q is a struct with the
%   fields {Q0, qdiag}, which define the process noise. See
%   C19FILT_KALMAN for their use.
%
%   On return, covZ is a covariance struct with fields:
%
%   Field         Dimensions
%   ------------------------------------------------------------------
%   lastPriorY    NSTATS-by-NNODES           Last prior state
%   lastPostY     NSTATS-by-NNODES           Last posterior state
%   lastPriorYCov (NSTATS-by-NNODES)^2       Last prior covariance
%   lastPostYCov  (NSTATS-by-NNODES)^2       Last posterior covariance
%
%   stdZ          size(G,1)-by-NTIME-by-NRATES
%   covZ          size(G,1)^2-by-NTIME-by-NRATES
%
%   The local state is understood as Y_loc = [I A E phi H W D R] and
%   the full state Y is just Y_loc replicated over NNODES. The output
%   observable is thus Z = G*Y, with covariance covZ = G*covY*G' and
%   stdZ = sqrt(diag(covZ)).
%
%   YDATA is KDATA-by-NNODES-by-MTIME (with MTIME <= NTIME) and
%   defines the measured signals (e.g., [H W D] where KDATA = 3) over
%   time. SLAB is a length MTIME vector defining what slab each data
%   point should be modeled as.
%
%   [Y,covY,resids] = C19FILT(...) additionally returns the residuals.
%
%   [Y,covY,resids,L] = C19FILT(...) additionally computes the marginal
%   (log-)likelihood L of the Kalman filter, see C19FILT_KALMAN. The
%   output matrix L is NTIME-by-NRATES.
%
%   [Y,covY,resids,L,S] = C19FILT(...) additionally returns the
%   innovation covariance S, where S(:,:,m,n) defines the covariance at
%   time index m and rate n.
%
%   See also INITC19FILT, C19FILT_KALMAN, GETC19SYST, GETC19OBS.

% H. Runvik 2021-02-18 (Adding Q input)
% R. Eriksson 2021-01-28 (merge with dynamical_beta branch)
% S. Engblom 2020-11-16 (Revision, output covY+G)
% S. Engblom 2020-10-31 (Revision, output L)
% S. Engblom 2020-10-22 (Revision, slab)
% H. Runvik 2020-10-20 (obsrates)
% H. Runvik, S. Engblom 2020-10-02 (H+P rather than H)
% S. Engblom 2020-09-18 (Major revision)
% H. Runvik 2020-09-18 (Revision, updated structure)
% S. Engblom 2020-09-17 (Revision, spatially connected model)
% S. Engblom 2020-09-02 (Revision, tspan added)
% H. Runvik 2020-08-19 (Extended number of states)
% S. Engblom 2020-04-29 (Revision)
% S. Engblom 2020-04-24 (Revision)
% S. Engblom 2020-04-15

% sizes
nstate = obsrates.nstate;
nlan = size(Ydata,2); % "nnodes"
nmeasured = numel(obsrates.states);
if isfield(obsrates,'indobs')
  nmeasured = nmeasured +numel(obsrates.indobs);
end

% read data & change into a more convenient format:
if nargin < 11
  % model error (Q0+qdiag*[...])
  Id = speye(nstate*nlan);
  Q.Q0 = Id;
  Q.qdiag = 0.05^2; % relative term, taken as +/- 5% error in rates
  % Alternatively one can increase the considired error in the
  % rates. We did however find the above to be optimal through trial
  % and error.
  % Q.qdiag = 0.01^2; % relative term, taken as +/- 1% error in rates
  % Q.qdiag = 0.1^2;% relative term, taken as +/- 10% error in rates
  % Q.qdiag = 0.25^2; % relative term, taken as +/- 25% error in rates
  if nargin < 10
    beta = [];
    if nargin < 9
      ntime = size(Ydata,3);
      if nargin < 8
        slab = ones(1,ntime);
        if nargin < 7
          error('Too few input arguments.');
        end
      end
    end
  end
end

if nargin < 12
  exception = [];
end
% effective transmission matrix
Deff = D-diag(sum(D,1));

% allocate output
nout = size(G,1);
if size(G,2) ~= nstate*nlan
  error('Output operator mismatch.');
end
nparams = size(rates.beta,2);
Z = zeros(nout,ntime,nparams);
covZ_aux = zeros(nout,nout,ntime);
covZ.stdZ = zeros(nout,ntime,nparams);
covZ.covZ = zeros(nout,nout,ntime,nparams); % not needed it seems
% (used to index out of (the posterior of) P, diag(COV)):
ixPdiag = find(repmat(eye(nout),[1 1 ntime]));
covZ.lastPriorY = zeros(nstate*nlan,nparams);
covZ.lastPostY = zeros(nstate*nlan,nparams);
covZ.lastPriorYCov = zeros(nstate*nlan,nstate*nlan,nparams);
covZ.lastPostYCov = zeros(nstate*nlan,nstate*nlan,nparams);

% output sizing
resids = zeros(size(Ydata,1)*size(Ydata,2),ntime,nparams);
if nargout > 3
  L = zeros(ntime,nparams);
  computeL = true;
  if nargout > 4
    S = zeros(nmeasured,nmeasured,ntime,nparams);
    computeS = true;
  else
    computeS = false;
  end
else
  computeL = false;
end

% logic to determine the effective slabs
slabloop = 1+[0 find(diff(slab) ~= 0) numel(slab)];
slabloop(end) = slabloop(end)-1; % (convention: simplifies code)
nloop = numel(slabloop)-1;
nslab = max(slab); % (note: must be contiguous)
F_loc = cell(1,nslab);
Qp = cell(1,nslab);
Qv = cell(1,nslab);
F = cell(1,nslab);

% switch to prediction by the end (optional)
mtime = size(Ydata,3);
if numel(slab) ~= mtime
  error('Data and slabs mismatch.');
end
if ntime > mtime
  % C19filt_kalman understands Nan's as missing data:
  Ydata = cat(3,Ydata,NaN(size(Ydata,1),size(Ydata,2),ntime-mtime));
  % extend the final slab correspondingly:
  slabloop(end) = slabloop(end)+ntime-mtime;
end

% network operator
T = zeros(nstate); T(4,4) = 1; % (variable phi)
Id_lan = speye(nlan);
Deff = kron(Deff,T);

% exception handling
if isempty(exception)
  exception.LB = -1e2;%-inf;
  exception.UB = 1e7;%inf;
  exception.SDFAC = 0.25;%0.5;
  exception.AbsMagn = 1e4;%1e6;
end


% loop over parameters
for n = 1:nparams
    % create all Kalman filters
    for m = 1:nslab
        % create new obsrates.indobspars for each sample of k_sew
        % used to construct matrix H below
        if isfield(rates, 'k_sew')
            obsrates.indobspars = {rates.k_sew(m,n)};
        end

        C19 = getC19syst(rates,m,n);
        KF = getC19filt(C19);
        F_loc{m} = KF.F;
        Qp{m} = KF.Qp;
        Qv{m} = KF.Qv;

    % "networkify" the model: this construct understands the state to be
    % ordered nstate states at a time, repeated nlan times on an
    % nstate*nlan long vector

    % network transmission model
    D = rates.lambda(m,n)*Deff;
    if any(size(F_loc{m}) ~= nstate)
      % annoying: cannot detect this until here
      error('Number of state variables mismatch.');
    end
    F{m} = kron(Id_lan,F_loc{m})+D;
    % note that the Qp has a bit of a special format with Q assembled on
    % the fly
  end

  % filter measurement models:

  % 1st phase (H,W,D), i.e., state (#5,#6,#7)
  if ~isfield(obsrates, 'R0_p2') % Equivalent to no k_sew amongst rates
    [H_loc,R0,rdiag,~,~] = getC19obs(obsrates);
    % measurement error model
    R0 = kron(Id_lan,R0);
    rdiag = repmat(rdiag,nlan,1);
  else
    [H_loc, R0, rdiag, R0_p2, rdiag_p2] = getC19obs(obsrates);
    % measurement error model
    R0 = kron(Id_lan,R0);
    rdiag = repmat(rdiag,nlan,1);
    R0_p2 = kron(Id_lan, R0_p2);
    rdiag_p2 = repmat(rdiag,nlan,1);
  end

  H = kron(Id_lan,H_loc);

  % initial state
  xx = zeros(nstate*nlan,1);
  PP = zeros(nstate*nlan);
  for k = 1:nlan
    yInit = Ydata(:,k,1);
    if isempty(Y0) && isempty(Y0Cov)
    [xx((1+(k-1)*nstate):k*nstate), ...
      PP((1+(k-1)*nstate):k*nstate,(1+(k-1)*nstate):k*nstate)] = ...
      initC19filt(F_loc{slab(slabloop(1))},H_loc,yInit,KF.CStates);
    else % use origo
       [xx((1+(k-1)*nstate):k*nstate), ...
      PP((1+(k-1)*nstate):k*nstate,(1+(k-1)*nstate):k*nstate)] = ...
      initC19filt(F_loc{slab(slabloop(1))},H_loc,yInit,KF.CStates,...
               Y0((1+(k-1)*nstate):k*nstate),...
               Y0Cov((1+(k-1)*nstate):k*nstate,(1+(k-1)*nstate):k*nstate));
    end
  end

  % assembly of constant part of Kalman filter
  KF.Q0 = Q.Q0;
  KF.qdiag = Q.qdiag;

  % exception handling
  KF.exception = exception;

  % measurements according to H
  KF.H = H;
  x = zeros(nstate*nlan,size(Ydata,3));
  x_prior = zeros(nstate*nlan,size(Ydata,3));
  P = zeros(nstate*nlan,nstate*nlan,size(Ydata,3));

  for j = 1:nloop
    slabix = slab(slabloop(j));
    KF.F = F{slabix};
    KF.Qp = Qp{slabix};
    KF.Qv = Qv{slabix};
      
    if exist('obsrates.sew_change', 'var') & j > obsrates.sew_change & exist('R0_p2','var')
      % second period constants:
      KF.R0 = R0_p2; 
      KF.rdiag = rdiag_p2;
    else
      % first period constants:
      KF.R0 = R0;
      KF.rdiag = rdiag;
    end

    dataslab = reshape(Ydata(:,:,slabloop(j):slabloop(j+1)),size(H,1),[]);
    if isempty(beta)
        betaslab=[];
    else
        betaslab.betaIdx=beta.betaIdx;
        betaslab.beta=beta.beta(slabloop(j):slabloop(j+1));
    end
    if computeL
      [x_,P_,resids(:,slabloop(j):slabloop(j+1),n), ...
       L(slabloop(j):slabloop(j+1),n),S_] = ...
          C19filt_kalman(KF,dataslab,xx,PP,betaslab);
      if computeS
        S(:,:,slabloop(j):slabloop(j+1),n)=S_;
      end
    else
      [x_,P_,resids(:,slabloop(j):slabloop(j+1),n)] = ...
          C19filt_kalman(KF,dataslab,xx,PP,betaslab);
    end

    % absorb posterior output (detail: by construction, the last entry
    % slabloop(j+1) is overwritten in the next round, except in the
    % final round)
    x(:,slabloop(j):slabloop(j+1)) = x_(:,:,2);
    x_prior(:,slabloop(j):slabloop(j+1)) = x_(:,:,1); % store prior 4 debug
    P(:,:,slabloop(j):slabloop(j+1)) = P_(:,:,:,2);

    % initial state for next slab: (note: *prior* is used)
    xx = x_(:,end,1);
    PP = P_(:,:,end,1);
  end

  % compute output states
  Z(:,:,n) = G*x;
  for i = 1:size(P,3)
    covZ_aux(:,:,i) = G*P(:,:,i)*G';
  end
  covZ.covZ(:,:,:,n) = covZ_aux;
  covZ.stdZ(:,:,n) = reshape(sqrt(covZ_aux(ixPdiag)),nout,ntime);

  % last prior/posterior:
  covZ.lastPriorY(:,n) = x_(:,end,1);
  covZ.lastPostY(:,n) = x_(:,end,2);
  covZ.lastPriorYCov(:,:,n) = P_(:,:,end,1);
  covZ.lastPostYCov(:,:,n) = P_(:,:,end,2);
end
