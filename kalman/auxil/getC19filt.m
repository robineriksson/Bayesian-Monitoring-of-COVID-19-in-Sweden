function KF = getC19filt(C19)
%GETC19FILT state- and covariance matrices for Covid-19 model.
%   KF = GETC19FILT(C19) produces a Kalman filter KF for the Covid-19
%   model C19, see GETC19SYST. KF is meant to be used by
%   C19FILT and C19FILT_KALMAN.
%
%   On return, KF contains four fields:
%
%   -F: State update matrix (i.e., x_{k+1} = F*x_k)
%
%   -Qp: State dependent covariance based on a Poisson assumption,
%   Qp(:,:,i) defines the contribution of the covariance from state
%   x_i according to Q_i = Qp(:,:,i)*x_i, Qpsum = sum(Q_i).
%   
%   -Qv: State dependent covariance assumed proportional to flows (i.e.
%   incidence) squared. Qv(i,k) defines the contribution from state k to
%   the variance of state i (no covariance between states included in this)
%   according to Q_i = sum_k(Qv(i,k)*x_k^2).
%
%   -CStates: logical array defining the cumulative states of the system
%   (states without outgoing transitions in the transition
%   matrix or additional cumulative states as defined by C19.CStates)
%
%   See also C19FILT, GETC19SYST, C19FILT_KALMAN.

% S. Engblom 2021-02-28 (Revision)
% H. Runvik 2020-11-18
% H. Runvik 2020-10-25
% H. Runvik 2020-10-21 (Revision)
% H. Runvik 2020-09-14

% convenience
transMat = C19.transMat;
dMat = C19.dMat;
CStates = C19.CStates;
DStates = C19.DStates;
AStates = C19.AStates;
IStates = C19.IStates;

% sizes
N = size(transMat,2); % "bulk"
Nc = numel(CStates);  % cumulative states
Ni = numel(IStates);  % incidense states
Ntot = N+Nc+Ni;

% checks
assert(all(size(dMat) == [numel(DStates) N]), ...
       'Incompatible sizes.');
assert(all(transMat(:,DStates) == 0,'all'), ...
       'Transition into deterministic state.');
assert(all(diag(transMat) == 0), ...
       'Transition into self.');
assert(all(transMat(AStates,:) == 0,'all'), ...
       'Transition from removed state.');
assert(all(dMat(:,AStates) == 0,'all'), ...
       'Transition from removed state into deterministic state.');
checkD = true(1,N); checkD(DStates) = false;
assert(all(checkD(CStates)) && all(checkD(IStates)), ...
       'Deterministic state does not support incidence/accumulation.');

% construct all cumulative states
allCStates = all(transMat == 0,2);
allCStates(AStates) = [];
KF.CStates = [allCStates; true(Nc,1); false(Ni,1)];

% start from no transitions
F = eye(N);
Qp = zeros(N,N,N);

% (1) bulk transitions
F = F+transMat'-diag(sum(transMat,2));
Qv = (transMat').^2+diag(sum(transMat.^2,2));
%Qv(DStates,DStates) = 0; % (not needed - done later)
for k1 = 1:N % FROM deterministic/stochastic state
  v = transMat(k1,:);
  Qp(:,:,k1) = diag(v);
  if checkD(k1) % FROM stochastic state
    Qp(k1,k1,k1) = Qp(k1,k1,k1)+sum(v,2); % add variance
    Qp(k1,:,k1) = Qp(k1,:,k1)-v;          % covariance
    Qp(:,k1,k1) = Qp(:,k1,k1)-v';         % covariance
  end
end

% (2) cumulative states
F = [F; F(CStates,:)];                % base state
F((N+1:N+Nc)+(CStates-1)*(N+Nc)) = 0; % remove coupling to base state
Qv = [Qv; Qv(CStates,:)];
Qv((N+1:N+Nc)+(CStates-1)*(N+Nc)) = 0;
Qp = [[Qp; zeros(Nc,N,N)] zeros(N+Nc,Nc,N)];
for k = 1:Nc
  k2 = 1:N;
  k2(CStates(k)) = [];
  Qp(N+k,:,k2) = Qp(CStates(k),:,k2);
  Qp(:,N+k,k2) = Qp(:,CStates(k),k2);
  Qp(N+k,N+k,k2) = Qp(CStates(k),CStates(k),k2);
end

% (3) incidence states
F = [F; F(IStates,:)];                 % base state
F((N+Nc+1:Ntot)+(IStates-1)*Ntot) = 0; % remove coupling to base state
Qv = [Qv; Qv(IStates,:)];
Qv((N+Nc+1:Ntot)+(IStates-1)*Ntot) = 0;
Qp = [[Qp; zeros(Ni,N+Nc,N)] zeros(N+Nc+Ni,Ni,N)];
for k = 1:Ni
  k2 = 1:N;
  k2(IStates(k)) = [];
  Qp(N+Nc+k,:,k2) = Qp(IStates(k),:,k2);
  Qp(:,N+Nc+k,k2) = Qp(:,IStates(k),k2);
  Qp(N+Nc+k,N+Nc+k,k2) = Qp(IStates(k),IStates(k),k2);
end

% (4) deterministic states
F(DStates,:) = dMat;
Qv(DStates,:) = dMat.^2;
Qv(DStates,DStates) = 0;

% (5) pad with zeros, remove unwanted states
F = [F zeros(Ntot,Nc+Ni)];
F((N+1:N+Nc)+(N:N+Nc-1)*Ntot) = 1; % "cumulative diagonal one's"
F(AStates,:) = [];
F(:,AStates) = [];

Na = numel(AStates); % "true" Ntot = Ntot-Na
Qp(AStates,:,:) = [];
Qp(:,AStates,:) = [];
Qp(:,:,AStates) = [];
Qp = cat(3,Qp,zeros(Ntot-Na,Ntot-Na,Nc+Ni));

Qv(AStates,:) = [];
Qv(:,AStates) = [];
Qv = [Qv zeros(Ntot-Na,Nc+Ni)];

% final assign
KF.F = F;
KF.Qp = Qp;
KF.Qv = Qv;
