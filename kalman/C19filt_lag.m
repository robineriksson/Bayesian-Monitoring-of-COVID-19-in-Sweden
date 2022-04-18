function [Z,covZ,ZLag,covZLag,L] = C19filt_lag(G,rates,D,obsrates,Ydata,slab,ntime,lag,Q)
%C19FILT_LAG (Lagged) Kalman filter for the EngEr Covid-19 model given data.
%   [Z,covZ,ZLag,covZLag,LAGLIKE] = C19FILT_LAG(...) returns the same
%   result [Z,covZ] as C19FILT, but with an additional prediction
%   [ZLag,covZLag] made LAG days ahead of each point in time, and the
%   associated (quasi-)log-likelihood LAGLIKE.
%
%   LAGLIKE is a cell of cells, 1 per SLAB in each rate combination in
%   RATES.
%
%   See also C19FILT.

% S. Engblom 2020-10-22 (Revision, slab)
% H. Runvik, S. Engblom 2020-10-02 (H+P rather than H)
% S. Engblom 2020-09-18 (Major revision)
% H. Runvik 2020-09-18 (Revision, updated structure)
% S. Engblom 2020-09-17 (Revision, spatially connected model)
% S. Engblom 2020-09-02 (Revision, tspan added)
% H. Runvik 2020-08-19 (Extended number of states)
% S. Engblom 2020-04-29 (Revision)
% S. Engblom 2020-04-24 (Revision)
% S. Engblom 2020-04-15

% read data & change into a more convenient format:
nlan = size(Ydata,2);
if nargin < 6, ntime = size(Ydata,3); slab=ones(1,ntime); lag=1; end
if nargin < 7, slab=ones(1,ntime); lag=1; end
if nargin < 8, lag=1; end
if nargin < 9
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
end

% effective transmission matrix
Deff = D-diag(sum(D,1));

% allocate output
nstate = obsrates.nstate;
nlan = size(Ydata,2); % "nnodes"
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
%ixPdiagLag = find(repmat(eye(nstate*nlan),[1 1 max(ntime,size(Ydata,3))]));
%YLag = zeros(nlan,nstate,max(ntime,size(Ydata,3)),nparams);
if nargout > 2
  if lag > 1
    computeLag = true;
    ZLag = zeros(nout,ntime,nparams);
    covZLag_aux = zeros(nout,nout,ntime);
    covZLag.stdZ = zeros(nout,ntime,nparams);
    covZLag.covZ = zeros(nout,nout,ntime,nparams);
  else
    computeLag = false;
    ZLag = NaN;
    covZLag = NaN;
  end

  if nargout > 4
    L = zeros(ntime,nparams);
    computeL = true;
  else
    computeL = false;
  end
end

% logic to determine the effective slabs
slabloop = 1+[0 find(diff(slab) ~= 0) numel(slab)];
slabloop(end) = slabloop(end)-1; % (convention: simplifies code)
nloop = numel(slabloop)-1;
nslab = max(slab); % (note: must be contiguous)
F_loc = cell(1,nslab);
Qp = cell(1,nslab);
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

% loop over parameters
for n = 1:nparams
  % create Kalman filters
  for m = 1:nslab
    C19 = getC19syst(rates,m,n);
    KF = getC19filt(C19);
    F_loc{m} = KF.F;
    Qp{m} = KF.Qp;

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

  % 1st phase (H+P,W,D), i.e., state (#5 + #8,#6,#7)
  [H_loc,R0,rdiag] = getC19obs(obsrates);
  H = kron(Id_lan,H_loc);

  % measurement error model
  R0 = kron(Id_lan,R0);
  rdiag = repmat(rdiag,nlan,1);

  % initial state
  xx = zeros(nstate*nlan,1);
  PP = zeros(nstate*nlan);
  for k = 1:nlan
    yInit = Ydata(:,k,1);
    [xx((1+(k-1)*nstate):k*nstate),...
     PP((1+(k-1)*nstate):k*nstate,(1+(k-1)*nstate):k*nstate)] = ...
        initC19filt(F_loc{slab(slabloop(1))},H_loc,yInit,KF.CStates,[],[]);
  end
  if computeLag
    xxLag = zeros(nlan*nstate,0);
    PPLag = zeros(nlan*nstate,nlan*nstate,0);
  end

  % base assembly
  KF.Q0 = Q.Q0;
  KF.qdiag = Q.qdiag;
  KF.R0 = R0;
  KF.rdiag = rdiag;

  % measurements according to H1
  KF.H = H;

  x = zeros(nstate*nlan,size(Ydata,3));
  P = zeros(nstate*nlan,nstate*nlan,size(Ydata,3));

  for j = 1:nloop
    slabix = slab(slabloop(j));
    KF.F = F{slabix};
    KF.Qp = Qp{slabix};

    dataslab = reshape(Ydata(:,:,slabloop(j):slabloop(j+1)),size(H,1),[]);
    if computeL
      [x_,P_,L(slabloop(j):slabloop(j+1),n)] = ...
        C19filt_kalman(KF,dataslab,xx,PP);
    else
      [x_,P_] = C19filt_kalman(KF,dataslab,xx,PP);
    end

    if computeLag
      %%%
      % Needs to consider the prediction ahead aswell for when mtime >
      % ntime. This is not really working as of yet, but should be fixed.
      %%%
      % make lag step forward predictions.
      if j < nloop
        tspanLag = 1:(size(dataslab,2));
        len = numel(tspanLag);
      else
        if ntime > size(Ydata,3)
          tspanLag = 1:(size(dataslab,2));
        else
          tspanLag = 1:(size(dataslab,2)-lag);
        end
        len = numel(tspanLag);
      end

      xLagFinal = zeros(nlan*nstate,len+lag);
      PLagFinal = zeros(nlan*nstate,nlan*nstate,len+lag);

      for i = tspanLag
        % look ahead lag days
        Ydata_empty = NaN(size(H,1),lag);
        [xLag,PLag] = C19filt_kalman(KF,Ydata_empty,x_(:,i,2),P_(:,:,i,2));
        % save final step in *Lag, as this is the lag steps prediction.
        ix=i+lag;
        xLagFinal(:,ix) = xLag(:,end,2);
        PLagFinal(:,:,ix) = PLag(:,:,end,2);
      end
      if j > 1
        xxLag = cat(2,xxLag, xLagFinal(:,(lag+2):end,:)); % not 100% certain about this.
        PPLag = cat(3,PPLag, PLagFinal(:,:,(lag+2):end,:));
      else
        xxLag = cat(2,xxLag, xLagFinal);
        PPLag = cat(3,PPLag, PLagFinal);
      end
    end

    % absorb posterior output (detail: by construction, the last entry
    % slabloop(j+1) is overwritten in the next round, except in the
    % final round)
    x(:,slabloop(j):slabloop(j+1)) = x_(:,:,2);
    P(:,:,slabloop(j):slabloop(j+1)) = P_(:,:,:,2);

    % initial state for next round:
    xx = x_(:,end,2);
    PP = P_(:,:,end,2);

  end


  % compute output states
  Z(:,:,n) = G*x;
  for i = 1:size(P,3)
    covZ_aux(:,:,i) = G*P(:,:,i)*G';
  end
  covZ.covZ(:,:,:,n) = covZ_aux;
  covZ.stdZ(:,:,n) = reshape(sqrt(covZ_aux(ixPdiag)),nout,ntime);
  if computeLag
    ZLag(:,:,n) = G*xxLag;
    for i = 1:size(P,3)
      covZLag_aux(:,:,i) = G*PPLag(:,:,i)*G';
    end
    covZLag.covZ(:,:,:,n) = covZ_aux;
    covZLag.stdZ(:,:,n) = reshape(sqrt(covZLag_aux(ixPdiag)),nout,ntime);
  end
end

% last prior/posterior:
covZ.lastPriorY = x_(:,end,1);
covZ.lastPostY = x_(:,end,2);
covZ.lastPriorYCov = P_(:,:,end,1);
covZ.lastPostYCov = P_(:,:,end,2);

if computeLag
  covZLag.lastPriorY = x_(:,end,1);
  covZLag.lastPostY = x_(:,end,2);
  covZLag.lastPriorYCov = P_(:,:,end,1);
  covZLag.lastPostYCov = P_(:,:,end,2);
end
